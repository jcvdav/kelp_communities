######################################################
#title#
######################################################
# 
# Purpose
#
######################################################

# Load packages
library(here)
library(vegan)
library(cowplot)
library(tidyverse)

kelp <- readRDS(here("data", "processed_data", "gebs_pisco_data.rds")) 
temps <- read_csv(here("data", "processed_data", "sst_by_location.csv"))


################################################################################

# BLatitudinal differences
abundance_matrix_latitude <- kelp %>% 
  filter(transect > 0,
         year == "Winter 2011 - 2012") %>% 
  group_by(source, location, site, level, transect, genus_species) %>% 
  summarize(abundance = sqrt(sum(abundance, na.rm = T))) %>% 
  ungroup() %>% 
  group_by(source, location, genus_species) %>% 
  summarize(abundance = mean(abundance, na.rm = T)) %>%
  ungroup() %>% 
  spread(genus_species, abundance, fill = 0)

metadata_abundance_matrix_latitude <- abundance_matrix_latitude %>% 
  select(source, location)

distance_matrix_latitude <- abundance_matrix_latitude %>% 
  select(-c(source, location)) %>% 
  vegdist(method = "bray")

set.seed(43)
latitudinal_nmds <- metaMDS(distance_matrix_latitude, trace = F)

latitudinal_nmds_data <- cbind(metadata_abundance_matrix_latitude, scores(latitudinal_nmds)) %>% 
  left_join(temps, by = "location")

latitudinal_nmds_plot <- ggplot(data = latitudinal_nmds_data, aes(x = NMDS1, y = NMDS2, fill = sstmean)) +
  geom_point(aes(shape = source), size = 4) +
  annotate(geom = "text", x = -0.1, y = -0.45, label = paste("2D Stress =", formatC(latitudinal_nmds$grstress, digits = 4, format = "f"))) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_gradientn(colors = colorRamps::matlab.like(100))+
  guides(fill = guide_colorbar(title = "Mean SST (ºC)", frame.colour = "black", ticks.colour = "black"),
         shape = guide_legend(title = "Program"))

latitudinal_nmds_plot2 <- ggplot(data = latitudinal_nmds_data, aes(x = NMDS1, y = NMDS2, color = sstmean)) +
  geom_text(aes(label = location)) +
  annotate(geom = "text", x = -0.1, y = -0.45, label = paste("2D Stress =", formatC(latitudinal_nmds$grstress, digits = 4, format = "f"))) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_gradientn(colors = colorRamps::matlab.like(100))+
  guides(color= guide_colorbar(title = "Mean SST (ºC)", frame.colour = "black", ticks.colour = "black"),
         shape = guide_legend(title = "Program"))




# Time differences
twice_sampled <- count(kelp, year, location) %>%
  distinct() %>%
  count(location) %>%
  filter(n > 1) %>% 
  pull(location) %>% 
  as.character()


# Time differences
abundance_matrix_time <- kelp %>% 
  filter(transect > 0,
         location %in% twice_sampled) %>% 
  group_by(source, year, location, site, level, transect, genus_species) %>% 
  summarize(abundance = sqrt(sum(abundance, na.rm = T))) %>% 
  ungroup() %>% 
  group_by(source, year, location, genus_species) %>% 
  summarize(abundance = mean(abundance, na.rm = T)) %>%
  ungroup() %>% 
  spread(genus_species, abundance, fill = 0)

metadata_abundance_matrix_time <- abundance_matrix_time %>% 
  select(source, year, location)

distance_matrix_time <- abundance_matrix_time %>% 
  select(-c(source, year, location)) %>% 
  vegdist(method = "bray")

set.seed(43)
nmds_time <- metaMDS(distance_matrix_time, trace = F)


mds_data <- cbind(metadata_abundance_matrix_time, scores(nmds_time))


time_nmds_plot <- ggplot(data = mds_data, aes(x = NMDS1, y = NMDS2, fill = year)) +
  geom_point(size = 4, aes(shape = source)) +
  geom_path(aes(group = location)) +
  annotate(geom = "text", x = -0.35, y = 0.35,
           label = paste("2D Stress =", formatC(nmds_time$grstress, digits = 4, format = "f"))) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("black", "white")) +
  guides(shape = "none",
         fill = guide_legend("Season", override.aes = list(shape = 21)))



leg_lat <- get_legend(latitudinal_nmds_plot)
leg_time <- get_legend(time_nmds_plot)

latitudinal_nmds_plot_sl <- latitudinal_nmds_plot +
  theme(legend.position = "none")
time_nmds_plot_sl <- time_nmds_plot +
  theme(legend.position = "none")

leg <- plot_grid(leg_lat, leg_time, ncol = 1)

nmds_plot <- plot_grid(latitudinal_nmds_plot_sl, time_nmds_plot_sl, ncol = 1, labels = "AUTO")

final_figure <- plot_grid(nmds_plot, leg, rel_widths = c(3.5, 1))

ggsave(plot = final_figure, filename = here("results", "img", "nmds_figure.png"),
       width = 6, height = 9)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################









################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################




# Simper analysis to test for community structure (2011-2012 data only)
simper_data <- kelp %>% 
  filter(transect > 0,
         year == "W 2011 - 2012") %>% 
  group_by(location, site, level, transect, genus_species) %>%	
  summarize(n = sqrt(sum(abundance))) %>%	
  ungroup() %>% 	
  spread(genus_species, n, fill = 0)	


comm <- simper_data %>% 	
  select_if(is.numeric) %>% 
  select(-transect) %>% 
  as.matrix()

sim <- simper(comm = comm,	
              group = simper_data$location,	
              permutations = 99,	
              parallel = 6)

res <- sim %>% 
  summary() %>% 
  map_dfr(magrittr::extract, .id = "pair") %>%
  rownames_to_column() %>%
  mutate(rowname = str_remove_all(rowname, "[:digit:]|[\\.+]")) %>% 
  rename(genus_species = rowname) %>% 
  as_tibble() %>% 
  group_by(pair) %>% 
  mutate(average_norm = average / sum(average)) %>% 
  ungroup() %>% 
  mutate(genus_species = fct_reorder(genus_species, average_norm, max))


simper_plot <- res %>% 
  filter(cumsum < 0.70, p < 0.05) %>% 
  ggplot(aes(x = genus_species, y = average_norm)) +	
  geom_boxplot(outlier.size = 0) +
  geom_jitter(size = 0.1, width = 0.25, height = 0) +
  theme_bw() +	
  theme(axis.text.y = element_text(face = "italic")) +
  coord_flip() +	
  labs(x = "Species", y = "Percent disimilarity") +	
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent)

ggsave(plot = simper_plot,
       filename = here("results", "img", "simper_plot.png"),
       width = 6, height = 9)



































