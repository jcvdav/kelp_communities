

library(here)
library(janitor)
library(readxl)
library(tidyverse)

gebs_data <- read.csv(here("data", "raw_data", "Peces_KelpForest_2011-2013.csv"),
                 stringsAsFactors = F) %>% 
  transform(location = reorder(location, latitude)) %>% 
  select(id, year, location, latitude, longitude, site, level, transect, kelp_density, genus_species, total_length, abundance) %>% 
  mutate(source = "GEBS") %>% 
  mutate(genus_species = case_when(genus_species == "Engraulix mordax/A. compressa" ~ "E. mordax/A. compressa",
                                   T ~ genus_species))

pisco_data <- read_xls(path = here("data", "raw_data", "PISCO_kelpforest_fish_ARV_211019.xls"), sheet = "DB2_Final", na = "NA") %>% 
  select(id, year,  location, latitude, longitude, site, level, transect, kelp_density, genus_species, total_length, abundance) %>% 
  mutate(location = case_when(location == "POINT_DUME" ~ "PDU",             
                              location == "LITTLE_IRISH_CEN" ~ "LIC",       
                              location == "SMI_HARRIS_PT_RESERVE_E" ~ "SHP",
                              location == "CANNERY_DC" ~ "CAN",             
                              location == "ANACAPA_MIDDLE_ISLE_CEN" ~ "ANA",
                              location == "SCI_CAVERN_POINT_E" ~ "SCI",     
                              location == "NAPLES_CEN" ~ "NAP")) %>% 
  mutate(source = "PISCO") %>% 
  mutate(genus_species = case_when(genus_species == "Engraulis mordax" ~ "E. mordax/A. compressa",
                                   genus_species == "Anchoa compressa" ~ "E. mordax/A. compressa",
                                   genus_species == "Sebastes serranoides,flavidus" ~ "Sebastes serranoides",
                                   genus_species %in% c("Sebastes YOY", "Stebastes complex YOY pteropodus") ~ "Sebastes spp",
                                   genus_species == "Anisotremus davidsoni" ~ "Anisotremus davidsonii",
                                   genus_species == "Halicoeres semicinctus" ~ "Halichoeres semicinctus",
                                   genus_species == "Sebastes dalli" ~ "Sebastes dallii",
                                   T ~ genus_species)) %>% 
  filter(!genus_species %in% c("ND", "Bathymasteridae", "Cottidae", "Sebastes spp"))


combined <- rbind(pisco_data, gebs_data) %>% 
  mutate(year = case_when(year <= 2012 ~ "Winter 2011 - 2012",
                          year == 2013 ~ "Winter 2013 - 2014")) %>% 
  mutate(location = fct_reorder(location, latitude)) %>% 
  group_by(location) %>% 
  mutate(latitude = mean(latitude),
         longitude = mean(longitude))

write_csv(x = combined, file = here("data", "processed_data", "gebs_pisco_data.csv"))
saveRDS(object = combined, file = here("data", "processed_data", "gebs_pisco_data.rds"))
