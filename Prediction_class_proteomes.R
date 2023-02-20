###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure XX -- ####
library(tidyverse)

## Load relevant data
load("./output/d.R")
load("./output/meta_pg.R")

## Calculate expression per class
read_csv("./input/kmeans_5_cluster.csv") %>%
  dplyr::select(-Slide) %>%
  mutate(cluster = ifelse(cluster == 5, 2,
                   ifelse(cluster == 4, 3,
                   ifelse(cluster == 3, 4,
                   ifelse(cluster == 2, 5, 1))))) -> NN_classes #Rename classes in spatially correct order


SA_incl_NN <- length(unique(NN_classes$cell_ID)) 

d %>%
  filter(cell_ID %in% NN_classes$cell_ID) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  mutate(completeness = n()/SA_incl_NN) %>%
  filter(completeness  == 1) %>%
  distinct(Protein) %>%
  pull(Protein) -> proteome_complete_NN

d %>%
  filter(Protein %in% proteome_complete_NN) %>%
  dplyr::select(Protein, cell_ID, int_core) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Protein) %>%
  left_join(NN_classes) %>%
  mutate(int_core = 2^int_core) %>%
  mutate(int_core = replace_na(int_core, 0)) %>%
  group_by(cluster, Protein) %>%
  summarise(median = mean(int_core, na.rm = T)) %>%
  left_join(meta_pg) -> d_classes

save(d_classes, file = "./Output/d_classes.R")