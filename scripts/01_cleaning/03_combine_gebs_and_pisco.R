

library(here)
library(janitor)
library(readxl)
library(tidyverse)

gebs_data <- read.csv(here("data", "raw_data", "Peces_KelpForest_2011-2013.csv"),
                 stringsAsFactors = F) %>% 
  transform(location = reorder(location, latitude)) %>% 
  select(id, year, location, latitude, longitude, site, level, transect, kelp_density, genus_species, total_length, abundance) %>% 
  mutate(source = "GEBS")

pisco_data <- read_xls(path = here("data", "raw_data", "PISCO_kelpforest_fish_ARV_211019.xls"), sheet = "DB2_Final", na = "NA") %>% 
  select(id, year,  location, latitude, longitude, site, level, transect, kelp_density, genus_species, total_length, abundance) %>% 
  mutate(location = case_when(location == "POINT_DUME" ~ "PDU",             
                              location == "LITTLE_IRISH_CEN" ~ "LIC",       
                              location == "SMI_HARRIS_PT_RESERVE_E" ~ "SMI",
                              location == "CANNERY_DC" ~ "CAN",             
                              location == "ANACAPA_MIDDLE_ISLE_CEN" ~ "ANA",
                              location == "SCI_CAVERN_POINT_E" ~ "SCI",     
                              location == "NAPLES_CEN" ~ "NAP")) %>% 
  mutate(source = "PISCO")

    


combined <- rbind(pisco_data, gebs_data) %>% 
  mutate(year = case_when(year <= 2012 ~ "Winter 2011 - 2012",
                          year == 2013 ~ "Winter 2013 - 2014"))

write_csv(x = combined, file = here("data", "processed_data", "gebs_pisco_data.csv"))
