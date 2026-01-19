# Load libraries
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(stringr)
library(tidyr)
library(ggplot2)

# Load and clean occurrence data
occ <- read_csv("data/com_data.csv")
analysis_df <- read_csv("output/analysis_df.csv")

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

# Create a grid-level presence table
grid_presence <- analysis_df %>%
  select(grid_id, source, presence) %>%
  distinct() %>%
  pivot_wider(
    names_from  = source,
    values_from = presence,
    values_fill = 0
  )

# Create the source-combination category
grid_presence <- grid_presence %>%
  mutate(
    source_combo = case_when(
      CS == 1 & Museum == 0 & Facebook == 0 ~ "Citizen science",
      CS == 0 & Museum == 1 & Facebook == 0 ~ "Museum",
      CS == 0 & Museum == 0 & Facebook == 1 ~ "Facebook",
      CS == 1 & Museum == 1 & Facebook == 0 ~ "Citizen science + Museum",
      CS == 1 & Museum == 0 & Facebook == 1 ~ "Citizen science + Facebook",
      CS == 0 & Museum == 1 & Facebook == 1 ~ "Museum + Facebook",
      CS == 1 & Museum == 1 & Facebook == 1 ~ "Citizen science + Museum + Facebook",
      TRUE ~ NA_character_
    )
  )

table(grid_presence$source_combo)

# Citizen science: 3812
# Facebook: 732
# Museum: 51
# Citizen science + Facebook: 436
# Citizen science + Museum: 17
# Museum + Facebook: 7
# Citizen science + Museum + Facebook: 9

# Join categories back to the spatial grid
grid_map <- grid %>%
  left_join(grid_presence, by = "grid_id")

grid_map_plot <- grid_map %>%
  filter(!is.na(source_combo))

# Get polygon
countries <- ne_countries(
  scale = "large",
  country = c("Bangladesh", "India", "Myanmar"),
  returnclass = "sf"
)

countries_utm <- st_transform(countries, crs = st_crs(crs_proj))

# Define a clear, colourblind-safe palette
source_cols <- c(
  "Citizen science" = "darkgoldenrod1",
  "Museum" = "tomato",
  "Facebook" = "blue",
  "Citizen science + Museum" = "#66C2A5",
  "Citizen science + Facebook" = "#8DA0CB",
  "Museum + Facebook" = "skyblue",
  "Citizen science + Museum + Facebook" = "black"
) 

# Plot
ggplot(grid_map_plot) +
  geom_sf(aes(fill = source_combo), colour = NA) +
  geom_sf(data = countries_utm, fill = NA, color = "grey90", size = 0.1) +  # outline on top
  scale_fill_manual(
    values = source_cols,
    na.value = "grey90",
    name = "Data source combination"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

ggsave("output/figures/spatial_bias_dist.png")

###############################################
# Prepare data in long format (lat + lon together)
grid_coords_long <- grid_coords %>%
  pivot_longer(
    cols = c(latitude, longitude),
    names_to = "coordinate",
    values_to = "value"
  ) %>%
  mutate(
    coordinate = recode(
      coordinate,
      latitude = "Latitude",
      longitude = "Longitude"
    )
  )


# Set the required order explicitly
grid_coords_long <- grid_coords_long %>%
  mutate(
    source_combo = factor(
      source_combo,
      levels = c(
        "Citizen science",
        "Museum",
        "Facebook",
        "Citizen science + Museum",
        "Citizen science + Facebook",
        "Museum + Facebook",
        "Citizen science + Museum + Facebook"
      )
    )
  )

# Combined latitudinal and longitudinal plot
ggplot(
  grid_coords_long,
  aes(
    x = source_combo,
    y = value,
    fill = source_combo
  )
) +
  geom_boxplot(
    outlier.alpha = 0.25,
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ coordinate,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = source_cols,
    guide = "none"
  ) +
  theme_bw() +
  theme(
    ## remove panel box
    panel.border = element_blank(),
    
    ## keep horizontal grid lines only
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    ## add axis lines back
    axis.line.x = element_line(colour = "black", linewidth = 0.5),
    axis.line.y = element_line(colour = "black", linewidth = 0.5),
    
    ## x-axis text
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    ## facet labels
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Data source combination",
    y = NULL
  )

ggsave("output/figures/spatial_bias_lat_lon.png")
