###########################
#### scDVP Figure Code ####
###########################

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Load relevant data
load("../output/variables/d.R")
load("../output/variables/meta_pg.R")

## Calculate expression per class
read_csv("../data/imaging/scDVP_kmeans5_clusters.csv") %>%
  dplyr::select(-Slide) -> NN_classes

SA_incl_NN <- length(unique(d$cell_ID)) 

d %>%
  dplyr::select(Protein, cell_ID, int_core) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Protein) %>%
  left_join(NN_classes) %>%
  mutate(int_core = 2^int_core) %>%
  mutate(int_core = replace_na(int_core, 0)) %>%
  group_by(cluster, Protein) %>%
  summarise(median = mean(int_core, na.rm = T)) %>%
  left_join(meta_pg) -> d_classes

save(d_classes, file = "../output/variables/d_classes.R")

