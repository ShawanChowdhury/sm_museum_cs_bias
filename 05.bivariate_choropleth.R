# Load libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(janitor)
library(biscale)
library(cowplot)

############################################################
# Import and clean occurrence data
############################################################
occ <- read_csv("data/com_data.csv", show_col_types = FALSE) %>%
  clean_names()

# Check your species column name.
# Change this if your file uses another name, e.g. "scientific_name"
species_col <- "species"

if (!species_col %in% names(occ)) {
  stop(
    paste0(
      "The species column '", species_col, "' was not found. ",
      "Please check the column names in occ and update species_col."
    )
  )
}

############################################################
# Set projection and create Bangladesh + India outline
############################################################
# Projected CRS in metres.
# This follows your original workflow.
crs_proj <- 32646  # UTM zone 46N

# Get Bangladesh + India polygons and dissolve into one outer polygon
countries_raw <- ne_countries(
  scale = "large",
  country = c("Bangladesh", "India"),
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_transform(crs_proj)

countries_utm <- countries_raw %>%
  st_union() %>%
  st_sf(geometry = .) %>%
  st_make_valid()

############################################################
# Clean and standardise occurrence data
############################################################
occ_clean <- occ %>%
  mutate(
    source = case_when(
      source == "HUMAN_OBSERVATION" ~ "CS",
      source == "PRESERVED_SPECIMEN" ~ "Museum",
      source == "Facebook" ~ "Facebook",
      TRUE ~ source
    )
  ) %>%
  filter(source %in% c("CS", "Museum", "Facebook")) %>%
  filter(
    !is.na(decimal_longitude),
    !is.na(decimal_latitude),
    !is.na(.data[[species_col]])
  )

# Convert occurrence records to sf points
occ_utm <- st_as_sf(
  occ_clean,
  coords = c("decimal_longitude", "decimal_latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(crs_proj)

# Keep only records inside Bangladesh + India
occ_utm <- occ_utm[lengths(st_intersects(occ_utm, countries_utm)) > 0, ]

############################################################
# Create 10 km grid and clip to Bangladesh + India outline
############################################################
grid <- st_make_grid(
  st_as_sfc(st_bbox(countries_utm)),
  cellsize = 10000,
  square = TRUE
) %>%
  st_sf(grid_id = seq_along(.))

# Keep only cells that intersect Bangladesh + India
grid <- grid[lengths(st_intersects(grid, countries_utm)) > 0, ]

# Clip grid to the outer country polygon
# This prevents grid cells from extending beyond the coastline/border
grid <- suppressWarnings(
  st_intersection(grid, countries_utm)
) %>%
  mutate(grid_id = row_number())

############################################################
# Assign occurrence records to grid cells
############################################################
occ_grid <- st_join(
  occ_utm,
  grid %>% select(grid_id),
  left = FALSE
)

############################################################
# Summarise richness and sampling effort
# Overall + separately for each data source
############################################################
# Function to create grid-level richness and effort summaries
summarise_grid_metrics <- function(data, label) {
  
  data %>%
    st_drop_geometry() %>%
    group_by(grid_id) %>%
    summarise(
      n_records = n(),
      n_species = n_distinct(.data[[species_col]]),
      .groups = "drop"
    ) %>%
    mutate(data_source = label)
}

# Overall map: all records combined
overall_metrics <- summarise_grid_metrics(
  data = occ_grid,
  label = "Overall"
)

# Source-specific maps
source_metrics <- occ_grid %>%
  st_drop_geometry() %>%
  group_by(grid_id, source) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(.data[[species_col]]),
    .groups = "drop"
  ) %>%
  mutate(
    data_source = recode(
      source,
      CS = "Citizen science",
      Museum = "Museum",
      Facebook = "Facebook"
    )
  ) %>%
  select(-source)

# Combine overall and source-specific summaries
all_metrics <- bind_rows(
  overall_metrics,
  source_metrics
)

############################################################
# Join summaries to grid and create manual bivariate classes
############################################################
grid_bivar <- grid %>%
  left_join(all_metrics, by = "grid_id") %>%
  filter(!is.na(data_source)) %>%
  mutate(
    richness_class = case_when(
      n_species == 1 ~ 1,
      n_species == 2 ~ 2,
      n_species >= 3 ~ 3,
      TRUE ~ NA_real_
    ),
    effort_class = case_when(
      n_records == 1 ~ 1,
      n_records %in% 2:3 ~ 2,
      n_records >= 4 ~ 3,
      TRUE ~ NA_real_
    ),
    bi_class = paste0(richness_class, "-", effort_class)
  )

# Check number of grid cells per map layer
print(table(grid_bivar$data_source))

# Citizen science        Facebook          Museum         Overall 
# 4125             794              46            4539 

# Check bivariate class distribution
print(table(grid_bivar$data_source, grid_bivar$bi_class))

# 1-1  1-2  1-3  2-2  2-3  3-2  3-3
# Citizen science 1951  243   26  581  135  162 1027
# Facebook         598   27    1  108   10   16   34
# Museum            22    5    8    4    0    1    6
# Overall         2201  248   25  658  130  177 1100

############################################################
# Function for bivariate choropleth map
############################################################
make_bivar_map <- function(data, title, outfile = NULL) {
  
  bb <- st_bbox(countries_utm)
  
  p_map <- ggplot() +
    geom_sf(
      data = data,
      aes(fill = bi_class),
      colour = NA,
      show.legend = FALSE
    ) +
    geom_sf(
      data = countries_utm,
      fill = NA,
      colour = "grey30",
      linewidth = 0.35
    ) +
    bi_scale_fill(
      pal = "DkBlue",
      dim = 3,
      na.value = "grey90"
    ) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      title = title,
      subtitle = "Richness: 1, 2, >=3 species; effort: 1, 2–3, >=4 records",
      x = NULL,
      y = NULL
    )
  
  p_legend <- bi_legend(
    pal = "DkBlue",
    dim = 3,
    xlab = "Higher observed richness",
    ylab = "Higher sampling effort",
    size = 8
  )
  
  p_final <- ggdraw() +
    draw_plot(p_map, 0, 0, 1, 1) +
    draw_plot(p_legend, 0.55, 0.70, 0.22, 0.22)
  
  if (!is.null(outfile)) {
    ggsave(
      filename = outfile,
      plot = p_final,
      width = 7,
      height = 8,
      dpi = 300
    )
  }
  
  return(p_final)
}

############################################################
# Overall bivariate choropleth map
############################################################
p_overall <- grid_bivar %>%
  filter(data_source == "Overall") %>%
  make_bivar_map(
    title = "Overall species richness and sampling effort",
    outfile = "output/figures/bivar_overall_richness_effort.png"
  )

p_overall

############################################################
# Source-specific bivariate choropleth maps
############################################################
source_levels <- c(
  "Citizen science",
  "Museum",
  "Facebook"
)

plot_list <- list()

for (src in source_levels) {
  
  dat_src <- grid_bivar %>%
    filter(data_source == src)
  
  if (nrow(dat_src) == 0) {
    message("Skipping ", src, ": no grid cells.")
    next
  }
  
  file_name <- src %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_$", "")
  
  p_src <- make_bivar_map(
    data = dat_src,
    title = paste0("Species richness and sampling effort: ", src),
    outfile = paste0("output/figures/bivar_", file_name, ".png")
  )
  
  plot_list[[src]] <- p_src
}

############################################################
# Save grid-level bivariate summary
############################################################
grid_bivar_summary <- grid_bivar %>%
  st_drop_geometry() %>%
  select(
    grid_id,
    data_source,
    n_records,
    n_species,
    richness_class,
    effort_class,
    bi_class
  )

write_csv(
  grid_bivar_summary,
  "output/grid_bivariate_summary_by_source.csv"
)
