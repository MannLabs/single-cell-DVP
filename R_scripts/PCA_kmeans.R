###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d_complete_heps.R")
load("../output/Variables/statTable_all.R")
load("../output/Variables/meta_distances.R")
# load("../output/Variables/meta_pg.R")
# source("./Functions/scale_dplyr.R")

## K-means clustering
set.seed(123)

clusters_heps <- as.data.frame(kmeans(t(d_complete_heps), centers = 5, iter.max = 1000, nstart = 50)[["cluster"]]) %>%
  rownames_to_column("cell_ID")

colnames(clusters_heps)[2] = "cluster"

statTable_all %>%
  distinct(bio_ID, label, cell_ID, run_ID, heps) %>%
  filter(cell_ID %in% colnames(d_complete_heps)) %>%
  arrange(cell_ID) %>%
  left_join(meta_distances) %>%
  mutate(run_number = as.numeric(str_replace(run_ID, ".*_", ""))) %>%
  drop_na(ratio) %>%
  left_join(clusters_heps) %>%
  dplyr::rename(Slide = bio_ID) %>%
  dplyr::select(Slide, Index, cell_ID, cluster) -> meta_heps_cluster

meta_heps_cluster %>%
  column_to_rownames("cell_ID") -> meta_heps_cluster_pca

p_cluster <- PCAtools::pca(d_complete_heps[,rownames(meta_heps_cluster_pca)], metadata = meta_heps_cluster_pca, removeVar = 0.1)

PCAtools::biplot(p_cluster ,
                 colby = 'cluster',
                 colLegendTitle = 'Cluster',
                 encircle = TRUE,
                 encircleFill = TRUE,
                 hline = 0, vline = c(-25, 0, 25),
                 legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0,
                 showLoadings = F, lab = NA)+
  scale_colour_viridis() -> plot_pca_kmeans

## Save Figure, table and R object
ggsave(file = "../output/Figures/PCA_kmeans.pdf", width = 7, height = 6, plot_pca_kmeans)
save(meta_heps_cluster, file = '../output/Variables/meta_heps_cluster.R')
write_tsv(meta_heps_cluster, file = "../output/Tables/kmeans_5clusters.csv")