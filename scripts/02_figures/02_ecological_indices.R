library(here)
library(cowplot)
library(tidyverse)

kelp <- readRDS(here("data", "processed_data", "gebs_pisco_data.rds"))

S_2013 <- kelp %>%
  group_by(year, source, location) %>%
  summarize(s2 = n_distinct(genus_species),
            s = n_distinct(genus_species[transect != 0])) %>%
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

D_2013 <- kelp %>% 
  filter(transect > 0) %>% 
  filter(!genus_species == "Engraulis mordax") %>% 
  group_by(year, source, location, site, level, transect) %>% 
  summarize(n = sum(abundance, na.rm = T)/60) %>% 
  group_by(year, source, location) %>% 
  summarize(D = mean(n, na.rm = T),
            sd = sd(n, na.rm = T)) %>% 
  ungroup() %>% 
  select(year, source, location, D, sd) %>% 
  mutate(year = as.character(year)) %>% 
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


L_2013 <- kelp %>% 
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
  mutate(year = as.character(year)) %>% 
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

legend <- cowplot::get_legend(L_2013)

L_2013 <- L_2013 +
  theme(legend.position = "none")

(index_2011_2013 <- plot_grid(plot_grid(S_2013, D_2013, L_2013, ncol = 3, labels = "AUTO"),
                         legend, ncol = 1, rel_heights = c(10, 1)))

ggsave(plot = index_2011_2013,
       filename = here("results", "img", "latitudinal_indices_2011_2013.png"),
       width = 8,
       height = 6)
