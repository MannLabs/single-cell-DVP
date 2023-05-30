###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 4A -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/img_fluovalues.R")
load("../output/variables/p_bins.R")
load("../output/variables/meta_distances_bins.R")
load("../output/variables/SA_incl_all.R")

# meta_binuc <- read_csv("../data/meta_binucleation.csv") %>%
#   filter(is.na(Notes) & (Classification == "Binuc" | Classification == "Mono"))  %>%
#   mutate(cell_ID = str_replace(cell_ID, "Dimethyl-n-", "target"))

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

ggsave(plot_FACS, file = "../Output/Figures/FACS.pdf", width = 5, height = 5)

# - Staining intensity by proteome cluster
p_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(img_fluovalues) %>%
  filter(Channel == "Alexa647" | Channel == "Alexa568") %>%
  group_by(Channel, bio_ID) %>%
  mutate(mean_scale = scale_df(Mean)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(abs(bin-9)), y = mean_scale, fill = Channel)) +
  geom_hline(yintercept = 0, lty = "dotted")+
  geom_boxplot()+
  scale_fill_manual(values = rev(viridis(4, option = "magma")[c(3,4)]))+
  theme_bw() +
  labs(x = "Proteome group", y = "Scaled fluorescence intensity") -> plot_fluo_by_proteome_bin

ggsave(plot_fluo_by_proteome_bin, file = "../Output/Figures/Proteome-bin_vs_fluorescence.pdf", width = 6, height = 5)

# - Staining intensity by spatial cluster
meta_distances_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(img_fluovalues) %>%
  filter(Channel == "Alexa647" | Channel == "Alexa568") %>%
  group_by(Channel, bio_ID) %>%
  mutate(mean_scale = scale_df(Mean)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(abs(bin-9)), y = mean_scale, fill = Channel)) +
  geom_hline(yintercept = 0, lty = "dotted")+
  geom_boxplot()+
  scale_fill_manual(values = rev(viridis(4, option = "magma")[c(3,4)]))+
  theme_bw() +
  labs(x = "Proteome group", y = "Scaled fluorescence intensity") -> plot_fluo_by_spatial_bin

ggsave(plot_fluo_by_spatial_bin, file = "../Output/Figures/Spatial-bin_vs_fluorescence.pdf", width = 6, height = 5)
# 
# # - Size of hepatocytes as histogram
# img_fluovalues %>%
#   left_join(meta_binuc) %>%
#   filter(Channel == "Alexa647") %>%
#   filter(cell_ID %in% SA_incl_all) -> d_histogram
# 
# ggplot()+ 
#   theme_classic()+
#   scale_x_continuous(limits = c(0,1500), breaks = seq(0,1500, by = 200), guide = guide_axis(n.dodge = 2))+
#   geom_density(data = d_histogram %>% drop_na(Classification), aes(x = Area, fill = Classification), alpha = 0.2) +
#   geom_density(data = d_histogram, aes(x = Area), color = "black", lty = "dotted")+
#   scale_fill_manual(values = viridis(4)[2:3])+
#   geom_vline(data = d_histogram %>% filter(Classification == "Binuc"), aes(xintercept = mean(Area)), size = 2, lty = "dashed", color = "#31688EFF")+
#   geom_vline(data = d_histogram %>% filter(Classification == "Mono"), aes(xintercept = mean(Area)), size = 2, lty = "dashed", color = "#35B779FF")+
#   geom_vline(data = d_histogram, aes(xintercept = mean(Area)), size = 1, lty = "dotted", color = "grey10")+
#   # add text labels for the mean values
#   annotate("text", x = mean(d_histogram %>% filter(Classification == "Binuc") %>% pull(Area)), y = 0.002, 
#            label = round(mean(d_histogram %>% filter(Classification == "Binuc") %>% pull(Area))), color = "#31688EFF", 
#            hjust = -0.2, size = 5)+
#   annotate("text", x = mean(d_histogram %>% filter(Classification == "Mono") %>% pull(Area)), y = 0.0023, 
#          label = round(mean(d_histogram %>% filter(Classification == "Mono") %>% pull(Area))), color = "#35B779FF", 
#          hjust = -0.2, size = 5)+
#   annotate("text", x = mean(d_histogram %>% pull(Area)), y = 0.0025, 
#            label = round(mean(d_histogram %>% pull(Area))), color = "grey10", 
#            hjust = -0.2, size = 5)+
#   labs(x= "Area (um2)") -> plot_size_distribution
# 
# ggplot(data = d_histogram, aes(x = Classification, y = Area, fill = Classification))+
#   geom_boxplot()+
#   theme_classic()+
#   scale_fill_manual(values = viridis(4)[2:3])+
#   scale_y_continuous(breaks = seq(0,1500, by = 200)) -> plot_size_distribution_boxplot
#   
# ggsave(plot_size_distribution, file = "../Output/Figures/Hepatocyte-size-distribution.pdf", width = 6, height = 5)
# ggsave(plot_size_distribution_boxplot, file = "../Output/Figures/Hepatocyte-size-distribution_boxplot.pdf", width = 6, height = 5)
#              