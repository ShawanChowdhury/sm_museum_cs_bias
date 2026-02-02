# Load libraries
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(janitor)

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
  filter(countryCode %in% c("IN", "BD"))

write_csv(gbif_data, "data/cleanedRecords_GBIF.csv")

# Cleaning memory
rm(gbif_data)

##############################
# Combining datasets
# Import datasets
gbif_data <- read_csv("data/cleanedRecords_GBIF.csv")
fb_data <- read_csv("data/dist_data_fb.csv")

gbif_data$order <- NULL
colnames(gbif_data)[6] <- "order"

# Facebook data is from 2013, so we are removing GBIF data before 2012
gbif_data <- gbif_data %>% 
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"))

# Subsetting by columns
gbif_data <- gbif_data %>% 
  dplyr::select(order, family, species, decimalLongitude, decimalLatitude, day, month, year, basisOfRecord)

fb_data <- fb_data %>% 
  dplyr::select(order, family, species, decimalLongitude, decimalLatitude, day, month, year)

# Adding data source
fb_data$source <- "Facebook"
colnames(gbif_data)[9] <- "source"

combined_data <- rbind(fb_data, gbif_data)

# Facebook data is from 2013, so we are removing GBIF data before 2012
combined_data <- combined_data %>% 
  filter(year > 2012 & year < 2025)

# Combining with species traits and IUCN status
trait_iucn <- read_csv("data/trait.csv")

trait_iucn_com <- dplyr::left_join(combined_data, trait_iucn, by = "species")

# Exporting data
write_csv(trait_iucn_com, "data/com_data.csv")

# Cleaning memory
rm(list = ls())

####################################
# Creating trait data table
# Load + rename trait columns
occ_raw <- read_csv("data/com_data.csv", show_col_types = FALSE)

colnames(occ_raw)[10:18] <- c("mass", "zonation", "activity", "foraging",
                              "venom", "n_gard", "range_size", "iucn_status",
                              "thrt_status")
# Filter data frame
occ <- occ_raw %>%
  clean_names() %>%
  select(-n_gard, -iucn_status)

# Standardise sources
occ <- occ %>%
  filter(source %in% c("HUMAN_OBSERVATION", "Facebook", "PRESERVED_SPECIMEN")) %>%
  mutate(source = case_when(
    source == "HUMAN_OBSERVATION"   ~ "CS",
    source == "PRESERVED_SPECIMEN" ~ "Museum",
    source == "Facebook"           ~ "Facebook",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(source))

# Make sure mass/range_size are numeric (safe if they’re already numeric)
occ <- occ %>%
  mutate(
    mass = suppressWarnings(as.numeric(mass)),
    range_size = suppressWarnings(as.numeric(range_size))
  )

# Convert to sf + filter to region
crs_proj <- "EPSG:32646"  # UTM 46N (metres)

occ_sf <- st_as_sf(occ, coords = c("decimal_longitude", "decimal_latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs_proj)

countries <- ne_countries(
  scale = "large",
  country = c("Bangladesh", "India"),
  returnclass = "sf"
) %>%
  st_transform(crs_proj) %>%
  st_union() %>%
  st_sf()

occ_sf <- occ_sf[st_within(occ_sf, countries, sparse = FALSE), ]

# Create 10 × 10 km grid over the region
grid <- st_make_grid(
  st_as_sfc(st_bbox(countries)),
  cellsize = 10000,  # 10 km in metres
  square = TRUE
) %>%
  st_sf(grid_id = seq_along(.))

# Keep only grid cells that intersect the region
grid <- grid[lengths(st_intersects(grid, countries)) > 0, ] %>%
  mutate(grid_id = row_number())

# Assign records to grids
occ_grid <- st_join(occ_sf, grid, join = st_within) %>%
  filter(!is.na(grid_id))

# Species → trait lookup (traits are species-level)
traits <- occ %>%
  distinct(species, mass, range_size, zonation, activity, foraging, venom, thrt_status)

# Build grid × source × species (species-weighted)
grid_species <- occ_grid %>%
  st_drop_geometry() %>%
  select(grid_id, source, species) %>%
  filter(!is.na(species)) %>%
  distinct(grid_id, source, species) %>%     # <- each species counts once per grid × source
  left_join(traits, by = "species")

# Diagnostics / intensity measures
n_records <- occ_grid %>%
  st_drop_geometry() %>%
  count(grid_id, source, name = "n_records")

n_species <- grid_species %>%
  count(grid_id, source, name = "n_species")

# Environmental covariates per grid
built <- rast("data/built_areas_crop.tif")
names(built) <- "built_areas"
hfp   <- rast("data/hfp_crop.tif")
names(hfp)   <- "hfp"
clim  <- rast("data/climate.tif")           # keeps whatever names are in that raster

# Merge rasters
vars <- c(built, hfp, clim)

# Change crs
grid_vect <- vect(grid) %>% project(crs(vars))

# Extract grid ID based information
grid_covariates <- terra::extract(
  vars,
  grid_vect,
  fun = mean,
  na.rm = TRUE
) %>%
  as_tibble() %>%
  mutate(grid_id = grid$grid_id) %>%
  select(-ID)

# Export output
write_csv(grid_covariates, "data/grid_covariates.csv")

# Trait summaries per grid × source (means across species within each grid)
# Continuous traits: mean across SPECIES in each grid × source
cont_vars <- grid_species %>%
  group_by(grid_id, source) %>%
  summarise(
    mean_mass       = mean(mass, na.rm = TRUE),
    mean_range_size = mean(range_size, na.rm = TRUE),
    .groups = "drop"
  )

# Categorical traits: proportion of SPECIES in each category within grid × source
activity_props <- grid_species %>%
  filter(!is.na(activity)) %>%
  count(grid_id, source, activity) %>%
  group_by(grid_id, source) %>%
  mutate(prop = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = activity, values_from = prop, values_fill = 0, names_prefix = "activity_")

zonation_props <- grid_species %>%
  filter(!is.na(zonation)) %>%
  count(grid_id, source, zonation) %>%
  group_by(grid_id, source) %>%
  mutate(prop = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = zonation, values_from = prop, values_fill = 0, names_prefix = "zonation_")

foraging_props <- grid_species %>%
  filter(!is.na(foraging)) %>%
  count(grid_id, source, foraging) %>%
  group_by(grid_id, source) %>%
  mutate(prop = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = foraging, values_from = prop, values_fill = 0, names_prefix = "foraging_")

venom_props <- grid_species %>%
  filter(!is.na(venom)) %>%
  count(grid_id, source, venom) %>%
  group_by(grid_id, source) %>%
  mutate(prop = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = venom, values_from = prop, values_fill = 0, names_prefix = "venom_")

threat_props <- grid_species %>%
  filter(!is.na(thrt_status)) %>%
  count(grid_id, source, thrt_status) %>%
  group_by(grid_id, source) %>%
  mutate(prop = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = thrt_status, values_from = prop, values_fill = 0, names_prefix = "threat_")

trait_grid_summary <- n_species %>%
  left_join(n_records,       by = c("grid_id", "source")) %>%
  left_join(cont_vars,       by = c("grid_id", "source")) %>%
  left_join(activity_props,  by = c("grid_id", "source")) %>%
  left_join(zonation_props,  by = c("grid_id", "source")) %>%
  left_join(foraging_props,  by = c("grid_id", "source")) %>%
  left_join(venom_props,     by = c("grid_id", "source")) %>%
  left_join(threat_props,    by = c("grid_id", "source"))

# Final analysis table (grid × source, fill zeros for missing source)
valid_grids <- occ_grid %>%
  st_drop_geometry() %>%
  distinct(grid_id)

analysis_df <- expand_grid(
  grid_id = valid_grids$grid_id,
  source  = c("CS", "Facebook", "Museum")
) %>%
  left_join(trait_grid_summary, by = c("grid_id", "source")) %>%
  left_join(grid_covariates,    by = "grid_id") %>%
  mutate(
    n_species = replace_na(n_species, 0L),
    n_records = replace_na(n_records, 0L)
  )

# Export data
write_csv(grid_species,        "output/grid_species_speciesweighted.csv")
write_csv(trait_grid_summary,  "output/trait_grid_summary_speciesbased.csv")
write_csv(analysis_df,         "output/analysis_df_speciesbased.csv")
