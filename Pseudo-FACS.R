###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure XX -- ####
library(tidyverse)

## Read relevant data
load("./output/img_fluovalues.R")
load("./output/p_bins.R")

## Additional function
scale_df<- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

## Plotting functions

# - Pseudo-FACS
img_fluovalues %>%
  dplyr::select(Mean, Channel, cell_ID) %>%
  spread(Channel, Mean) %>%
  #summarise(AF647Pos = Alexa647 > 500, AF568Pos = Alexa568 > 275, n = Alexa568 > 0) %>%
  ggplot(aes(x = Alexa568, y = Alexa647))+
  geom_point(color = "darkred", alpha = 0.5, size = 2)+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10() -> plot_FACS

ggsave(plot_FACS, file = "./Output/Figures/FACS.pdf", width = 5, height = 5)

# - Staining intensity by proteome cluster
p_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(img_fluovalues) %>%
  filter(Channel == "Alexa647" | Channel == "Alexa568") %>%
  group_by(Channel, bio_ID) %>%
  mutate(mean_scale = scale_df(Mean)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(bin), y = mean_scale, fill = Channel)) +
  geom_hline(yintercept = 0, lty = "dotted")+
  geom_boxplot()+
  scale_fill_manual(values = rev(viridis(4, option = "magma")[c(3,4)]))+
  theme_bw() +
  labs(x = "Proteome group", y = "Scaled fluorescence intensity") -> plot_fluo_by_bin

ggsave(plot_fluo_by_bin, file = "./Output/Figures/Proteome-bin_vs_fluorescence.pdf", width = 6, height = 5)