
library(here)
library(raster)
library(rnaturalearth)
library(tidyverse)

kelp <- read.csv(here("data", "processed_data", "gebs_pisco_data.csv"),
                 stringsAsFactors = F) %>% 
  transform(location = reorder(location, latitude)) %>% 
  filter(!genus_species == "ND") %>% 
  mutate(year = ifelse(year == 2012, 2011, year),
         year = factor(year))

pts <- select(kelp, location, latitude, longitude, source) %>% 
  distinct()

mex_coast <- ne_countries(country = c("Mexico", "United States of America"), returnclass = "sf", scale = 50)

e <- extent(-126, -106, 24,  40)

# SST data from: http://gmed.auckland.ac.nz/download.html
sst <- raster(here("data", "raw_data", "sstmean", "sstmean.asc")) %>% 
  crop(e) %>% 
  as.data.frame(xy = T)

map <- ggplot() +
  geom_raster(data = sst, aes(x = x, y = y, fill = layer), alpha = 0.9) +
  geom_contour(data = sst, aes(x = x, y = y, z = layer), color = "black", binwidth = 1) +
  geom_sf(data = mex_coast) + 
  geom_point(data = pts, aes(x = longitude, y = latitude, shape = source), size = 2) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(26, 38)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-124, -108)) +
  scale_fill_gradientn(colors = colorRamps::matlab.like(100))+
  guides(fill = guide_colorbar("Mean SST",
                               frame.colour = "black",
                               ticks.colour = "black"),
         shape = guide_legend("Monitoring\nprogram")) +
  labs(x = "", y = "")

ggsave(plot = map, filename = here("results", "img", "map.png"),
       width = 6, height = 5)
