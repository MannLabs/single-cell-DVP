###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure XX -- ####
library(tidyverse)

## Read relevant data
load("./output/d_all_norm.R")

## Plotting function

d_all_norm %>%
  filter(Protein %in% c("P10854", "P27661", "P68433", "P62806")) -> d_histones

ggplot(data = d_histones, aes(x = Symbol, y = int_core, color = Symbol, shape = heps)) +
  geom_jitter()+
  theme_bw() +
    scale_color_manual(values = viridis(6)[2:5]) +
    scale_shape_manual(values = c(21,19))-> plot_histones

ggsave(plot_histones, file = "./Output/Figures/Histone_intensity.pdf", width = 6, height = 5)