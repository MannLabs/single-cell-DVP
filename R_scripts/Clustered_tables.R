###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 3A -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

## -- Halpern et al, RNASeq data, 9 spatial clusters
## Define number of classes
classes = 9

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> distance_bins_tmp

distance_bins_tmp %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(distance_bins_tmp) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> distance_bins

d %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Protein) %>%
  mutate(int_core = ifelse(is.na(int_core), 0, int_core)) %>%
  left_join(distance_bins %>% rownames_to_column("cell_ID") %>% drop_na(bin)) %>%
  left_join(meta_pg) %>%
  group_by(Protein, ENSEMBL, Symbol, bin) %>%
  drop_na(bin) %>%
  summarise(int = log2(mean(2^int_core))) %>%
  drop_na(bin) -> table_proteome_to_RNASeq

write_tsv(as.data.frame(table_proteome_to_RNASeq), "../output/Tables/Proteome_to_RNASeq_9spatialbins.tsv")

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/p_bins.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

## -- Ben-Moshe et al., FASC/Proteomics data, 8 PCA clusters
## Define number of classes
classes = 9

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% rownames(p_bins)) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = rownames(p_bins), pc1 = p_bins$pc1) %>%
  mutate(range = cut_interval(pc1, n = classes))  -> p_bins_tmp

p_bins_tmp %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(p_bins_tmp) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> p_bins

d %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Protein) %>%
  mutate(int_core = ifelse(is.na(int_core), 0, int_core)) %>%
  left_join(meta_pg) %>%
  left_join(p_bins %>% rownames_to_column("cell_ID") %>% drop_na(bin)) %>%
  group_by(Protein, Symbol, ENSEMBL, bin) %>%
  summarise(int = log2(mean(2^int_core))) %>%
  drop_na(bin) -> table_proteome_to_FACS

write_tsv(as.data.frame(table_proteome_to_FACS), "../output/Tables/Proteome_to_FACS_9PCbins.tsv")