###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/")
load("../output/variables/")

d_heatmap_binned_sorted %>%
  rownames_to_column("Protein") %>%
  gather(cluster_anchorprotein, int_mean, !Protein) %>%
  mutate(int_mean = 2^int_mean)
left_join(meta_pg) -> d_heatmap_binned_meta

write_csv(d_heatmap_binned_meta, "./RNAseq/221122_scDVP_Protoemics-to-RNAseq__distance_cluster_9.csv")

d_heatmap %>%
  rownames_to_column("Protein") %>%
  gather(cell_ID, int_core, !Protein) %>%
  drop_na(int_core) %>%
  left_join(meta_distances_bin) %>%
  drop_na(bin) %>%
  group_by(Protein, bin) %>%
  summarise(n = n()) %>%
  filter(n > 5) %>%
  left_join(d_heatmap %>%
              rownames_to_column("Protein") %>%
              gather(cell_ID, int_core, !Protein) %>%
              #left_join(p_bins %>% rownames_to_column("cell_ID")) %>%
              left_join(meta_distances_bin) %>%
              drop_na(bin) %>%
              group_by(Protein, bin) %>%
              summarise(median = median(int_core, na.rm = T))) %>%
  left_join(meta_pg) %>%
  drop_na(Symbol) %>%
  left_join(limma_8 %>% dplyr::select(Symbol, logFC)) %>%
  mutate(zonated = ifelse(logFC > 0, "CV", "PV")) %>%
  mutate(int = 2^median) %>%
  group_by(Symbol, bin, zonated) %>%
  summarise(median = median(int), sd = sd(int)) %>%
  group_by(Symbol) %>%
  mutate(sum = sum(median)) %>%
  mutate(ratio = median/sum) %>%
  dplyr::select(Symbol, ratio, bin) %>%
  dplyr::rename(cluster_anchorprotein = bin) -> d_heatmap_binned_ratio

write_csv(d_heatmap_binned_ratio, "./RNAseq/221122_scDVP_Protoemics-to-RNAseq__distance_cluster_9__RATIOS.csv")