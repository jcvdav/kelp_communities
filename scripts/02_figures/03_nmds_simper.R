
library(here)
library(vegan)
library(cowplot)
library(tidyverse)


# All sites Winter 2011 - 2012
abundance_matrix_2011 <- kelp %>% 
  filter(transect > 0,
         year == "Winter 2011 - 2012") %>% 
  group_by(source, location, site, level, transect, genus_species) %>% 
  summarize(abundance = sqrt(sum(abundance, na.rm = T))) %>% 
  ungroup() %>% 
  group_by(source, location, genus_species) %>% 
  summarize(abundance = mean(abundance, na.rm = T)) %>%
  ungroup() %>% 
  spread(genus_species, abundance, fill = 0)

metadata_abundance_matrix_2011 <- abundance_matrix_2011 %>% 
  select(source, location)

distance_matrix_2011 <- abundance_matrix_2011 %>% 
  select(-c(source, location)) %>% 
  vegdist(method = "bray")

set.seed(43)
mds <- metaMDS(data_c_cedros_samples, trace = F)

stress <- paste("2D Stress =", formatC(mds$grstress, digits = 4, format = "f"))

(nmds_change <- cbind(metadata_abundance_matrix_2011, scores(mds)) %>% 
    ggplot(aes(x = NMDS1, y = NMDS2, fill = location)) +
    geom_point(aes(shape = source), size = 4) +
    annotate(geom = "text", x = -0.2, y = 0.35, label = stress) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    scale_shape_manual(values = c(21, 24)) +
    scale_fill_viridis_d() +
    guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)),
           shape = guide_legend(title = "Source", order = 0)))


################
simper_data <- kelp %>% 
  filter(transect > 0,
         year == "Winter 2011 - 2012") %>% 
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
  map_dfr(extract) %>%
  rownames_to_column() %>%
  mutate(rowname = str_remove_all(rowname, "[:digit:]")) %>% 
  rename(genus_species = rowname)

# cusums <- map_dfc(sim, "cusum")	%>% 
#   cbind(map_df(sim, "cusum")) %>% 	
#   gather(pair, cumsum, -species)	
# 
# overall <- map_df(sim, "overall") %>% 	
#   gather(pair, overall)	

res %>% 
  filter(cumsum < 0.95) %>% 
  gather(pair, average, -genus_species) %>%	
  # left_join(cusums, by = c("species", "pair")) %>%
  # left_join(overall, by = "pair") %>%
  # filter(cumsum <= 0.90) %>%
  transform(species = reorder(species, average)) %>% 	
  mutate(pair = str_replace(pair, "_", " vs. "),	
         average = average / overall) %>% 
  ggplot(aes(x = species, y = average, fill = average)) +	
  geom_violin() +	
  theme_bw() +	
  coord_flip() +	
  labs(x = "Spp", y = "Avg % disim") +	
  theme(legend.position = "none") +	
  scale_fill_viridis_d()
