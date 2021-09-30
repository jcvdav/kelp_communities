# Create table of presence / absence by site

library(here)
library(tidyverse)

kelp_presence_absence <- read.csv(here("data", "Peces_KelpForest_2011-2013.csv"),
                 stringsAsFactors = F) %>% 
  transform(location = reorder(location, latitude)) %>% 
  filter(!year == 2012) %>% 
  group_by(location, genus_species) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(n = ifelse(n > 0, 1, 0)) %>% 
  spread(genus_species, n, fill = 0)

write.csv(x = kelp_presence_absence,
          file = here("docs", "community", "tables", "kelp_presence_absence.csv"),
          row.names = F)
