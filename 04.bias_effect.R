# Load libraries
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyr)
library(broom)

# Import data
analysis_df <- read_csv("output/analysis_df.csv")

# Pivot to long format for plotting
id_cols <- c("grid_id", "source", "presence")

analysis_long <- analysis_df %>%
  pivot_longer(
    cols = -all_of(id_cols),
    names_to = "predictor",
    values_to = "value"
  )

# Export data
write_csv(analysis_long, "output/analysis_long.csv")

# Create the Overall rows
analysis_long_overall <- analysis_long %>%
  group_by(grid_id, predictor) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    presence = as.integer(any(presence == 1)),
    source = "Overall",
    .groups = "drop"
  )

# Bind back to the long dataframe
analysis_long2 <- bind_rows(
  analysis_long,
  analysis_long_overall
)

analysis_long2 <- analysis_long2 %>%
  mutate(
    source = factor(source, levels = c("Overall", "CS", "Facebook", "Museum"))
  )

# Fit linear models (one per predictor)
lm_results <- analysis_long2 %>%
  filter(!is.na(value)) %>%
  group_by(predictor) %>%
  do(
    tidy(
      lm(scale(value) ~ source, data = .),
      conf.int = TRUE
    )
  ) %>%
  ungroup() %>%
  filter(term != "(Intercept)")

# Convert model output into your plotting dataframe
bias_df <- lm_results %>%
  mutate(
    source = gsub("^source", "", term),
    bias = estimate,
    bias_low = conf.low,
    bias_high = conf.high,
    p_value = p.value,
    significance = ifelse(p.value < 0.05, "Significant", "Not significant")
  ) %>%
  select(
    predictor,
    source,
    bias,
    bias_low,
    bias_high,
    p_value,
    significance
  )

# Make significance a factor
bias_df <- bias_df %>%
  mutate(
    significance = factor(
      significance,
      levels = c("Not significant", "Significant")
    )
  )

# Predictor label cleaning
bias_df <- bias_df %>%
  mutate(
    predictor_clean = predictor %>%
      str_replace_all("_", " ") %>%
      str_replace_all("&", " & ")
  )

# Export dataframe
write_csv(bias_df, "output/bias_df_lm.csv")

# Create a helper dataframe with predictor positions
predictor_order <- c(
  "activity Cathemeral",
  "activity Crepuscular",
  "activity Diurnal",
  "activity Nocturnal",
  "foraging Active foraging",
  "foraging Mixed",
  "foraging Sit and Wait",
  "mean mass",
  "mean range size",
  "venom Yes",
  "venom No",
  "zonation Aquatic",
  "zonation Arboreal",
  "zonation Cryptic",
  "zonation Fossorial",
  "zonation Saxicolous",
  "zonation Semi-Aquatic",
  "zonation Terrestrial",
  "zonation Arboreal & Saxicolous",
  "zonation Arboreal & Terrestrial",
  "zonation Fossorial & Terrestrial",
  "zonation Arboreal & Saxicolous & Terrestrial",
  "threat Threatened",
  "threat Non-threatened",
  "threat Data Deficient",
  "threat Not Evaluated",
  "mean elev",
  "mean temp",
  "mean rain",
  "mean built areas",
  "mean hfp"
)

# Custom predictor order
bias_df <- bias_df %>%
  mutate(
    predictor_clean = factor(
      predictor_clean,
      levels = rev(predictor_order)
    )
  )


hline_df <- bias_df %>%
  distinct(predictor_clean) %>%
  mutate(
    y = as.numeric(predictor_clean)
  )

# Plot
ggplot(
  bias_df,
  aes(
    x = bias,
    y = predictor_clean,
    colour = source,
    shape = significance
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey90") +
  
  geom_hline(
    data = hline_df,
    aes(yintercept = y),
    inherit.aes = FALSE,
    colour = "grey90",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  geom_point(size = 2.8) +
  
  geom_errorbarh(
    aes(xmin = bias_low, xmax = bias_high),
    height = 0.25
  ) +
  
  scale_colour_manual(
    values = c(
      CS = "darkgoldenrod1",
      Facebook = "blue",
      Museum = "tomato"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Not significant" = 1,  # open circle
      "Significant"     = 16  # filled circle
    )
  ) +
  
  labs(
    x = "Bias (standardised deviation from the combined data)",
    y = NULL,
    colour = "Data source",
    shape = "Statistical significance"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 3),
    legend.position = "none"
  )


# Export output
ggsave("output/figures/bias.png")