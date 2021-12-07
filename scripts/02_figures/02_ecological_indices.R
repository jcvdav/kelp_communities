######################################################
#ecological_indices
######################################################
# 
# Calculate ecological indices 
#
######################################################

# SET UP #######################################################################
library(here)
library(cowplot)
library(tidyverse)

kelp <- readRDS(here("data", "processed_data", "gebs_pisco_data.rds"))

# PROCESSING ###################################################################

# I first calculate ecological indicators (put them in a table) and then
# plot them. Finally, I put them all together and export them.
# I also export the figures.


# Species richness 
S_data <- kelp %>%
  group_by(year, source, location) %>%
  summarize(s2 = n_distinct(genus_species),
            s = n_distinct(genus_species[transect != 0]))

S <- S_data %>% 
  ggplot(aes(x = location, y = s, fill = year)) +
  geom_pointrange(aes(ymin = s, ymax = s2, shape = source), position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.position = "none") +
  labs(x = "Site", y = "Richness (num. spp)") +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 24))


# Density
D_data <- kelp %>% 
  filter(transect > 0) %>% 
  filter(!genus_species == "Engraulis mordax") %>% 
  group_by(year, source, location, site, level, transect) %>% 
  summarize(n = sum(abundance, na.rm = T)/60) %>% 
  group_by(year, source, location) %>% 
  summarize(D = mean(n, na.rm = T),
            sd = sd(n, na.rm = T)) %>% 
  ungroup() %>% 
  select(year, source, location, D, sd) %>% 
  mutate(year = as.character(year))

D <- D_data %>% 
  ggplot(aes(x = location, y = D, fill = year)) +
  geom_pointrange(aes(ymin = D - sd, ymax = D + sd, shape = source), position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.position = "none") +
  labs(x = "", y = expression(Density~(log[10](Org~m^{-2})))) +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 24))


L_data <- kelp %>% 
  filter(transect > 0) %>% 
  group_by(year, source, location, site, level, transect, genus_species) %>% 
  summarize(ni = sum(abundance, na.rm = T)) %>% 
  filter(ni > 0) %>% 
  group_by(year, source, location, site, level, transect) %>% 
  mutate(N = sum(ni)) %>%  
  ungroup() %>% 
  group_by(year, source, location, site, level, transect) %>%
  summarize(Li = 1 - sum((ni * (ni - 1)) / (N * (N - 1)))) %>% 
  ungroup() %>% 
  group_by(year, source, location) %>% 
  summarize(L = mean(Li, na.rm = T),
            sd = sd(Li, na.rm = T)) %>% 
  ungroup() %>% 
  select(year, source, location, L, sd)%>% 
  mutate(year = as.character(year))

L <- L_data %>% 
  ggplot(aes(x = location, y = L, fill = year)) +
  geom_pointrange(aes(ymin = L - sd, ymax = L + sd, shape = source), position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent",
                                         size=0.5,
                                         linetype = "solid", 
                                         colour = "black")) +
  labs(x = "", y = "Simpson's index") +
  scale_fill_manual(values = c("black", "white")) +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend("Program"),
         fill = guide_legend("Season", override.aes = list(shape = 21)))

legend <- cowplot::get_legend(L)

L <- L +
  theme(legend.position = "none")

(index_2011_2013 <- plot_grid(plot_grid(S, D, L, ncol = 3, labels = "AUTO"),
                         legend, ncol = 1, rel_heights = c(10, 1)))

index_data <- S_data %>% 
  left_join(D_data, by = c("year", "source", "location")) %>% 
  left_join(L_data, by = c("year", "source", "location")) %>% 
  select(-contains("sd")) %>% 
  rename(richness_transct_0 = s2,
         standard_richness = s,
         Density = D,
         Simpson = L)

# EXPORTS ######################################################################

# Export figure
ggsave(plot = index_2011_2013,
       filename = here("results", "img", "latitudinal_indices_2011_2013.png"),
       width = 8,
       height = 6)

# Export table
write_csv(x = index_data,
          file = here("data", "processed_data", "ecological_indices_by_location_and_year.csv"))
