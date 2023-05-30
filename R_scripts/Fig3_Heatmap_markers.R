###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 3D -- ####

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
markers <- data.frame(Symbol = c("Ass1", "Asl", "Cps1", "Arg1", "Glul", "Cyp2e1", "Oat"),
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
  gather(Symbol, int_scale, 3:ncol(.)) %>%
  mutate(int_scale = ifelse(int_scale > 2, 2, ifelse(int_scale < -2, -2, int_scale)))-> d_heat

ggplot(d_heat) +
  geom_segment(aes(x = as.numeric(ratio), xend = as.numeric(ratio), y = 0, yend = 1, col = as.numeric(int_scale)), size = 1) +
  scale_color_viridis(option = "inferno",na.value = "grey80", ) +
  theme_classic()+
  facet_wrap(.~Symbol, ncol = 1) -> plot_markers

ggplot(data = d_heat,aes(x = ratio, y = int_scale))+
  geom_point()+
  facet_wrap(.~Symbol, ncol = 1)+
  geom_smooth() +
  geom_hline(yintercept = 0)+
  theme_bw() -> plot_scatter

ggsave(plot_markers, file = "../Output/Figures/Heatmap_markers.pdf", width = 6, height = 6)
