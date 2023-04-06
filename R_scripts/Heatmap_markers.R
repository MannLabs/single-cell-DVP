###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")
source("./Functions/scale_dplyr.R")

## Markers to consider
markers <- data.frame(Symbol = c("Ass1", "Asl", "Otc", "Cps1", "Arg1", "Glul", "Cyp2e1", "Oat"),
                      POI = "POI") %>%
  left_join(meta_pg)

d %>%
  left_join(markers) %>%
  drop_na(POI) %>%
  left_join(meta_distances) %>%
  drop_na(ratio) %>%
  dplyr::select(cell_ID, Symbol, int_core, ratio) %>%
  group_by(Symbol) %>%
  mutate(int_scale = scale_dplyr(int_core)) %>%
  ungroup() %>%
  dplyr::select(-int_core) %>%
  spread(Symbol, int_scale) %>%
  gather(Symbol, int_scale, 3:ncol(.)) -> d_heat

ggplot(d_heat) +
  geom_segment(aes(x = as.numeric(ratio), xend = as.numeric(ratio), y = 0, yend = 1, col = as.numeric(int_scale)), size = 1) +
  scale_color_viridis(option = "inferno",na.value = "grey80") +
  theme_classic()+
  facet_wrap(.~Symbol, ncol = 1) -> plot_markers

ggsave(plot_markers, file = "../Output/Figures/Heatmap_markers.pdf", width = 6, height = 6)
