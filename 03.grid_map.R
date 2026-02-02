# Load libraries
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(stringr)
library(tidyr)
library(ggplot2)
library(janitor)

# Load and clean occurrence data
occ <- read_csv("data/com_data.csv", show_col_types = FALSE) %>% clean_names()
analysis_df <- read_csv("output/analysis_df_speciesbased.csv", show_col_types = FALSE)

### Use a standard CRS (metres), avoids km/metre confusion
crs_proj <- 32646  # EPSG:32646 (UTM zone 46N, metres)

# Get polygon (Bangladesh + India only) and transform to UTM
countries_utm <- ne_countries(
  scale = "large",
  country = c("Bangladesh", "India"),
  returnclass = "sf"
) %>%
  st_transform(crs_proj) %>%
  st_union() %>%
  st_sf()

# Convert points to UTM (metres)
to_utm <- function(df) {
  st_as_sf(df, coords = c("decimal_longitude", "decimal_latitude"), crs = 4326) %>%
    st_transform(crs_proj)
}

cs_utm     <- occ %>% filter(source == "HUMAN_OBSERVATION")   %>% mutate(source = "CS")     %>% to_utm()
fb_utm     <- occ %>% filter(source == "Facebook")            %>%                           to_utm()
museum_utm <- occ %>% filter(source == "PRESERVED_SPECIMEN")  %>% mutate(source = "Museum") %>% to_utm()

inside_region <- function(x) x[st_within(x, countries_utm, sparse = FALSE), ]
cs_utm     <- inside_region(cs_utm)
fb_utm     <- inside_region(fb_utm)
museum_utm <- inside_region(museum_utm)

# Build the grid from the countries
grid <- st_make_grid(
  st_as_sfc(st_bbox(countries_utm)),
  cellsize = 10000,
  square = TRUE
) %>%
  st_sf(grid_id = seq_along(.))

# Filter grid to Bangladesh+India
grid <- grid[lengths(st_intersects(grid, countries_utm)) > 0, ] %>%
  mutate(grid_id = row_number())

# Presence distribution
grid_presence <- analysis_df %>%
  select(grid_id, source, n_species) %>%
  distinct() %>%
  mutate(presence = as.integer(n_species > 0)) %>%
  select(-n_species) %>%
  pivot_wider(names_from = source, values_from = presence, values_fill = 0)

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

# Total: 4539
# Citizen science: 4125
# Facebook: 794
# Museum: 46

# Citizen science: 3711 
# Facebook: 394 
# Museum: 18 
# Citizen science + Facebook: 388 
# Citizen science + Museum: 16 
# Museum + Facebook: 2 
# Citizen science + Museum + Facebook: 10

grid_map_plot <- grid %>%
  left_join(grid_presence, by = "grid_id") %>%
  filter(!is.na(source_combo))

# Specifying colours
source_cols <- c(
  "Citizen science" = "darkgoldenrod1",
  "Museum" = "tomato",
  "Facebook" = "blue",
  "Citizen science + Museum" = "#66C2A5",
  "Citizen science + Facebook" = "#8DA0CB",
  "Museum + Facebook" = "skyblue",
  "Citizen science + Museum + Facebook" = "black"
)

# Crop plotting window to Bangladesh+India bbox (removes excess blank space)
bb <- st_bbox(countries_utm)

ggplot() +
  geom_sf(data = countries_utm, fill = NA, colour = "grey80", linewidth = 0.2) +
  geom_sf(data = grid_map_plot, aes(fill = source_combo), colour = NA) +
  scale_fill_manual(values = source_cols, na.value = "grey90", name = "Data source combination") +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggsave("output/figures/spatial_bias_dist.png")

###############################################
# Centroids as POINT geometry (one per grid cell)
centroids_sf <- st_centroid(st_geometry(grid_map_plot))

# Make an sf object with those centroids + your attributes
centroids_tbl <- st_sf(
  grid_id      = grid_map_plot$grid_id,
  source_combo = grid_map_plot$source_combo,
  geometry     = centroids_sf,
  crs          = st_crs(grid_map_plot)
)

# Convert to lon/lat and extract coordinates
centroids_ll <- st_transform(centroids_tbl, 4326)
xy <- st_coordinates(centroids_ll)

grid_coords <- centroids_ll %>%
  st_drop_geometry() %>%
  mutate(
    longitude = xy[, 1],
    latitude  = xy[, 2]
  ) %>%
  select(grid_id, source_combo, latitude, longitude)

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
