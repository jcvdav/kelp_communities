

library(here)
library(janitor)
library(readxl)
library(tidyverse)

gebs_data <- read.csv(here("data", "raw_data", "Peces_KelpForest_2011-2013.csv"),
                 stringsAsFactors = F) %>% 
  transform(location = reorder(location, latitude)) %>% 
  filter(!year == 2012)

pisco_data <- read_xls(path = here("data", "raw_data", "PISCO_kelpforest_fish_ARV_211019.xls"), sheet = "DB1_Select-Filter")
