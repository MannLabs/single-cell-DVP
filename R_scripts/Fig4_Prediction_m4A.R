###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Load additional functions and environment
source("./Functions/fct_xml_to_polygon.R")
load("../output/variables/d_classes.R")

prediction_m4A <- read_csv("../output/Tables/shape_probability_m4A.csv") %>%
  dplyr::select(-1) %>%
  dplyr::rename(Index = `Shape Index`) %>%
  gather(cluster, probability, 1:5) %>%
  mutate(cluster = as.numeric(cluster) + 1) %>%
  right_join(d_classes) %>%
  mutate(median_weighted = probability * median) %>%
  dplyr::select(Index, median_weighted, Protein) %>%
  group_by(Index, Protein) %>%
  mutate(int_weighted = sum(median_weighted)) %>%
  mutate(bio_ID = "m4A")

library(XML)
library(sf)
library(RColorBrewer)
library(ggpubr)
library(platetools)

polygons_m4A <- xml_to_polygon("../data/imaging/EXP-221117_scDVP_m4A_PREDICTION_237.xml")

level_of <- data.frame(Index = c(1:length(polygons_m4A))) %>%
  left_join(prediction_m4A) %>%
  filter(Protein == "P16460") %>%
  distinct(Index, .keep_all = T)

st_sf(polygons_m4A, ID = as.factor(c(1:length(polygons_m4A))), level_of) %>%
  
  ggplot()+
  geom_sf(aes(fill = log2(int_weighted)))+
  scale_fill_viridis_c()+
  theme(legend.position="none") +
  theme_classic()+
  scale_y_continuous(limits = c(132000,145000))+
  scale_x_continuous(limits = c(0, -14000)) -> plot_stsf_Ass1

rm(level_of)

level_of <- data.frame(Index = c(1:length(polygons_m4A))) %>%
  left_join(prediction_m4A) %>%
  filter(Protein == "P15105") %>%
  distinct(Index, .keep_all = T)

st_sf(polygons_m4A, ID = as.factor(c(1:length(polygons_m4A))), level_of) %>%
  
  ggplot()+
  geom_sf(aes(fill = log2(int_weighted)))+
  scale_fill_viridis_c()+
  #geom_sf_text(aes(label = as.character(ID)), colour = "black") +
  theme(legend.position="none") +
  theme_classic()+
  scale_y_continuous(limits = c(132000,145000))+
  scale_x_continuous(limits = c(0, -14000)) -> plot_stsf_Glul

rm(level_of)

level_of <- data.frame(Index = c(1:length(polygons_m4A))) %>%
  left_join(prediction_m4A) %>%
  filter(Protein == "Q05421") %>%
  distinct(Index, .keep_all = T)

st_sf(polygons_m4A, ID = as.factor(c(1:length(polygons_m4A))), level_of) %>%
  
  ggplot()+
  geom_sf(aes(fill = log2(int_weighted)))+
  scale_fill_viridis_c()+
  theme(legend.position="none") +
  theme_classic()+
  scale_y_continuous(limits = c(132000,145000))+
  scale_x_continuous(limits = c(0, -14000)) -> plot_stsf_Cyp2e1

ggsave(plot_stsf_Cyp2e1, file = "../output/Figures/Spatialmap_Cyp2e1.pdf", width = 5, height = 5)
ggsave(plot_stsf_Ass1, file = "../output/Figures/Spatialmap_Ass1.pdf", width = 5, height = 5)
ggsave(plot_stsf_Glul, file = "../output/Figures/Spatialmap_Glul.pdf", width = 5, height = 5)