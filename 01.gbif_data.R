# Loading required packages
library(tidyverse)
library(rgbif)

# Species name
species <- read_csv("data/dist_data_fb.csv")
species <- unique(species$species)

# match the names 
gbif_taxon_keys <- species %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

# download the data
occ_download(
  pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"),
  format = "SIMPLE_CSV",
  user= "",pwd= "",email= ""
)

# Check download status
occ_download_wait('0080320-251120083545085')

# Citation
# GBIF Occurrence Download https://www.gbif.org/occurrence/download/0080323-251120083545085 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2026-01-07
