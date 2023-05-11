###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 2B -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../Output/Variables/d_all_norm.R")

## Plotting function

d_all_norm %>%
  filter(Protein %in% c("P10854", "P27661", "P68433", "P62806")) -> d_histones

d_histones %>%
  filter(heps == TRUE) %>%
  group_by(Symbol) %>%
  mutate(is.max = max(int_core) == int_core,is.min = min(int_core) == int_core) %>%
  mutate(ismaxmin = is.max + is.min) %>%
  right_join(d_histones) %>%
  mutate(mark = ifelse(heps == FALSE, "endo", ifelse(ismaxmin == 0, "null", "minmax"))) -> d_histones_to_plot

ggplot(d_histones_to_plot, aes(x = Symbol, y = int_core, color = Symbol, shape = mark)) +
  geom_point(position = position_jitter(seed = 1))+
  theme_bw() +
  scale_color_manual(values = viridis(6)[2:5]) +
  scale_shape_manual(values = c(21,17,19)) -> plot_histones

ggsave(plot_histones, file = "../Output/Figures/Figure-2_Histone_intensity.pdf", width = 6, height = 5)
