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
library(glmmTMB)
library(broom.mixed)

# Import data
analysis_df <- read_csv("output/analysis_df_speciesbased.csv")

# Pivot to long format for plotting
id_cols <- c("grid_id", "source", "n_species")

analysis_long <- analysis_df %>%
  select(-n_records) %>% 
  pivot_longer(
    cols = -all_of(id_cols),
    names_to = "predictor",
    values_to = "value"
  )

# Export data
write_csv(analysis_long, "output/analysis_long.csv")

analysis_long2 <- analysis_long %>%
  mutate(source = factor(source, levels = c("CS", "Facebook", "Museum"))) %>%
  ### NEW: drop undefined trait rows (n_species==0 -> trait NAs)
  filter(n_species > 0, !is.na(value)) %>%
  ### NEW: drop predictors with no variation (these cause NA coefficients/p-values)
  group_by(predictor) %>%
  filter(n_distinct(value) > 1) %>%
  ungroup() %>%
  ### NEW: stabilise proportion predictors using logit transform (still simple)
  mutate(
    is_prop = str_detect(predictor, "^(activity|zonation|foraging|venom|threat)_"),
    value_t = if_else(
      is_prop,
      qlogis(pmin(pmax(value, 0.001), 0.999)),  # clamp then logit
      value
    )
  )

# Fit linear models (one per predictor)
lm_results <- analysis_long2 %>%
  group_by(predictor) %>%
  do(
    tidy(
      lm(scale(value_t) ~ source, data = ., weights = n_species),
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
  select(predictor, source, bias, bias_low, bias_high, p_value, significance) %>%
  ### NEW: remove any remaining NAs so they don’t show in plot/legend
  filter(!is.na(bias), !is.na(significance), source %in% c("Facebook", "Museum")) %>%
  mutate(
    significance = factor(significance, levels = c("Not significant", "Significant")),
    predictor_clean = predictor %>%
      str_replace_all("_", " ") %>%
      str_replace_all("&", " & ")
  )

bias_df <- bias_df %>%
  mutate(
    predictor_clean = predictor %>%
      str_replace_all("_", " ") %>%
      str_replace_all("&", " & "),
    predictor_clean = case_when(                 # <-- add this mapping
      predictor_clean == "hfp"           ~ "mean hfp",
      predictor_clean == "built areas"   ~ "mean built areas",
      predictor_clean == "annual precip" ~ "mean rain",
      predictor_clean == "elevation"     ~ "mean elev",
      TRUE ~ predictor_clean
    )
  )

# Export dataframe
write_csv(bias_df, "output/bias_df_weighted_lm_vsCS.csv")

# # Make significance a factor
# bias_df <- bias_df %>%
#   mutate(significance = factor(significance, levels = c("Not significant", "Significant")))
# 
# 
# # Predictor label cleaning
# bias_df <- bias_df %>%
#   mutate(
#     predictor_clean = predictor %>%
#       str_replace_all("_", " ") %>%
#       str_replace_all("&", " & ")
#   )
# 
# write_csv(bias_df, "output/bias_df_lm.csv")

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
  filter(predictor_clean %in% predictor_order) %>%         # <-- keeps only the ones you order
  mutate(predictor_clean = factor(predictor_clean, levels = rev(predictor_order)))

hline_df <- bias_df %>%
  distinct(predictor_clean) %>%
  mutate(y = as.numeric(predictor_clean))

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
  # 0 is the “no difference vs CS” reference 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey90") +
  geom_point(size = 2.8) +
  geom_errorbarh(aes(xmin = bias_low, xmax = bias_high), height = 0.25) +
  scale_colour_manual(values = c(Facebook = "blue", Museum = "tomato")) +
  scale_shape_manual(values = c("Not significant" = 1, "Significant" = 16)) +
  labs(
    x = "Bias (standardised difference vs Citizen Science; weighted by n_species)",
    y = NULL,
    colour = "Data source",
    shape = "Statistical significance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) + 
  scale_y_discrete(expand = expansion(mult = c(0.03, 0.03))) + 
  geom_hline(data = hline_df, aes(yintercept = y - 0.5), inherit.aes = FALSE,
             linetype = "dashed", colour = "grey90", linewidth = 0.3) +
  geom_hline(data = hline_df, aes(yintercept = y + 0.5), inherit.aes = FALSE,
             linetype = "dashed", colour = "grey90", linewidth = 0.3)

# Export output
ggsave("output/figures/bias.png")
