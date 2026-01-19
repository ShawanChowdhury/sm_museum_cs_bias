# Load libraries
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)

##############################
# Cleaning GBIF data
# Reading GBIF data
gbif_data <- data.table::fread("data/gbif/gbif.csv")

# Removing blank cells
gbif_data <- gbif_data[!(is.na(gbif_data$species) | gbif_data$species == ""),]
gbif_data <- gbif_data[!(is.na(gbif_data$decimalLongitude) | gbif_data$decimalLongitude == ""),]
gbif_data <- gbif_data[!(is.na(gbif_data$decimalLatitude) | gbif_data$decimalLatitude == ""),]

# Removing duplicated records
gbif_data <- gbif_data[!duplicated(gbif_data),]

# Selecting relevant countries
gbif_data <- gbif_data %>% 
  filter(countryCode %in% c("IN", "BD", "MM"))

write_csv(gbif_data, "data/cleanedRecords_GBIF.csv")

# Cleaning memory
rm(gbif_data)

##############################
# Combining datasets
# Import datasets
gbif_data <- read_csv("data/cleanedRecords_GBIF.csv")
fb_data <- read_csv("data/dist_data_fb.csv")

# Facebook data is from 2013, so we are removing GBIF data before 2012
gbif_data <- gbif_data %>% 
  filter(year > 2012 & year < 2025 & basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"))

gbif_data$order <- NULL
colnames(gbif_data)[6] <- "order"

# Subsetting by columns
gbif_data <- gbif_data %>% 
  dplyr::select(order, family, species, decimalLongitude, decimalLatitude, day, month, year, basisOfRecord)

fb_data <- fb_data %>% 
  dplyr::select(order, family, species, decimalLongitude, decimalLatitude, day, month, year)

# Adding data source
fb_data$source <- "Facebook"
colnames(gbif_data)[9] <- "source"

combined_data <- rbind(fb_data, gbif_data)

# Combining with species traits and IUCN status
trait_iucn <- read_csv("data/trait.csv")

trait_iucn_com <- dplyr::left_join(combined_data, trait_iucn, by = "species")

# Exporting data
write_csv(trait_iucn_com, "data/com_data.csv")

# Cleaning memory
rm(list = ls())


####################################
# Creating trait data table
# Load and clean occurrence data
occ <- read_csv("data/com_data.csv")

# Define sources and create “overall” sampling frame
occ_all <- bind_rows(
  occ,
  occ %>% mutate(source = "Overall")
)

cs_data     <- occ_all %>% filter(source == "HUMAN_OBSERVATION") %>% mutate(source = "CS")
fb_data     <- occ_all %>% filter(source == "Facebook")
museum_data <- occ_all %>% filter(source == "PRESERVED_SPECIMEN") %>% mutate(source = "Museum")
ov_data     <- occ_all %>% filter(source == "Overall")

# Spatial conversion and projection
crs_proj <- "+proj=utm +zone=46 +datum=WGS84 +units=km +no_defs"

to_utm <- function(df) {
  st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    st_transform(crs_proj)
}

cs_utm     <- to_utm(cs_data)
fb_utm     <- to_utm(fb_data)
museum_utm <- to_utm(museum_data)
ov_utm     <- to_utm(ov_data)

# Define study region and filter points
countries <- ne_countries(
  scale = "large",
  country = c("Bangladesh", "India", "Myanmar"),
  returnclass = "sf"
) %>%
  st_transform(crs_proj) %>%
  st_union() %>%
  st_sf()

inside_region <- function(x) x[st_within(x, countries, sparse = FALSE), ]

cs_utm     <- inside_region(cs_utm)
fb_utm     <- inside_region(fb_utm)
museum_utm <- inside_region(museum_utm)
ov_utm     <- inside_region(ov_utm)

# Create 10 × 10 km grid (conditioned on overall data)
grid <- st_make_grid(
  st_as_sfc(st_bbox(ov_utm)),
  cellsize = c(10, 10),
  square = TRUE
) %>%
  st_sf(grid_id = seq_along(.))

# Assign points to grids
assign_grid <- function(points) {
  st_join(points, grid, join = st_within) %>%
    filter(!is.na(grid_id))
}

cs_grid     <- assign_grid(cs_utm)
fb_grid     <- assign_grid(fb_utm)
museum_grid <- assign_grid(museum_utm)
ov_grid     <- assign_grid(ov_utm)

# Binary grid × source presence (conditioning step)
valid_grids <- ov_grid %>%
  st_drop_geometry() %>%
  distinct(grid_id)

presence_df <- expand_grid(
  grid_id = valid_grids$grid_id,
  source  = c("CS", "Facebook", "Museum")
) %>%
  left_join(
    bind_rows(
      cs_grid     %>% st_drop_geometry() %>% distinct(grid_id) %>% mutate(source = "CS", presence = 1),
      fb_grid     %>% st_drop_geometry() %>% distinct(grid_id) %>% mutate(source = "Facebook", presence = 1),
      museum_grid %>% st_drop_geometry() %>% distinct(grid_id) %>% mutate(source = "Museum", presence = 1)
    ),
    by = c("grid_id", "source")
  ) %>%
  mutate(presence = replace_na(presence, 0))

# Species–trait lookup table
colnames(occ)[10:18] <- c("mass", "zonation", "activity", "foraging", 
                          "venom", "n_gard", "range_size", "iucn_status", 
                          "thrt_status")
occ$n_gard <- NULL
occ$iucn_status <- NULL

traits <- occ %>%
  distinct(
    species,
    mass,
    zonation,
    activity,
    foraging,
    venom,
    range_size,
    thrt_status
  ) %>%
  clean_names()

# Attach species to grid × source and summarise traits
grid_species <- bind_rows(
  cs_grid     %>% st_drop_geometry() %>% dplyr::select(grid_id, species) %>% mutate(source = "CS"),
  fb_grid     %>% st_drop_geometry() %>% dplyr::select(grid_id, species) %>% mutate(source = "Facebook"),
  museum_grid %>% st_drop_geometry() %>% dplyr::select(grid_id, species) %>% mutate(source = "Museum")
) %>%
  left_join(traits, by = "species")

# Extract grid-level environmental covariates
built <- rast("data/built_areas_crop.tif")
hfp   <- rast("data/hfp_crop.tif")
clim  <- rast("data/climate.tif")

vars <- c(built, hfp, clim)

grid_vect <- vect(grid) %>%
  project(crs(vars))

terraOptions(threads = parallel::detectCores() - 4)

grid_covariates <- terra::extract(
  vars,
  grid_vect,
  fun = mean,
  na.rm = TRUE
) %>%
  as_tibble() %>%
  mutate(grid_id = grid$grid_id) %>%
  select(-ID)

# Merge with the original dataset
grid_merge <- grid_species %>%
  left_join(grid_covariates, by = "grid_id")

# Export output
write_csv(grid_merge, "output/grid_merge.csv")

# Continuous variables
cont_vars <- grid_merge %>%
  group_by(grid_id, source) %>%
  summarise(
    mean_mass = mean(mass, na.rm = TRUE),
    mean_range_size = mean(range_size, na.rm = TRUE),
    mean_built_areas = mean(built_areas),
    mean_hfp = mean(hfp),
    mean_temp = mean(mean_temp),
    mean_rain = mean(annual_precip),
    mean_elev = mean(elevation),
    .groups = "drop"
  )

# Categorical traits (proportions)
cat_props <- function(df, var, prefix) {
  df %>%
    filter(!is.na(.data[[var]])) %>%
    count(grid_id, source, .data[[var]]) %>%
    group_by(grid_id, source) %>%
    mutate(prop = n / sum(n)) %>%
    select(-n) %>%
    pivot_wider(
      names_from = all_of(var),
      values_from = prop,
      values_fill = 0,
      names_prefix = prefix
    )
}

trait_grid_summary <- cont_vars %>%
  left_join(cat_props(grid_merge, "activity", "activity_"), by = c("grid_id", "source")) %>%
  left_join(cat_props(grid_merge, "zonation", "zonation_"), by = c("grid_id", "source")) %>%
  left_join(cat_props(grid_merge, "foraging", "foraging_"), by = c("grid_id", "source")) %>%
  left_join(cat_props(grid_merge, "venom", "venom_"), by = c("grid_id", "source")) %>%
  left_join(cat_props(grid_merge, "thrt_status", "threat_"), by = c("grid_id", "source"))

# Join presence and traits
analysis_df <- presence_df %>%
  left_join(trait_grid_summary, by = c("grid_id", "source"))

# Export output
write_csv(analysis_df, "output/analysis_df.csv")