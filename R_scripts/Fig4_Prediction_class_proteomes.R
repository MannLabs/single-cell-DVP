###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Load relevant data
load("../output/variables/d.R")
load("../output/variables/meta_pg.R")

## Calculate expression per class
read_csv("../output/Tables/kmeans_5clusters.csv") %>%
  dplyr::select(-Slide) -> NN_classes


SA_incl_NN <- length(unique(d$cell_ID)) 

d %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  filter(cell_ID %in% NN_classes$cell_ID) %>%
  drop_na(int_core) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  gather(cell_ID, int_core, !Protein) %>%
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

save(d_classes, file = "../output/variables/d_classes.R")

