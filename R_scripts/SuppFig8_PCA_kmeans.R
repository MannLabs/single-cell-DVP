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
load("../output/Variables/d.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/p_heps.R")
# # load("../output/Variables/meta_pg.R")
# # source("./Functions/scale_dplyr.R")
# 
# # library(kpodclustr)
# # 
# # d %>%
# #   filter(cell_ID %in% SA_incl_all) %>%
# #   dplyr::select(cell_ID, int_core, Protein) %>%
# #   spread(cell_ID, int_core) %>%
# #   column_to_rownames("Protein") -> d_wide
# # 
# # kpod <- kpod(as.matrix(t(d_wide)), k = 5)
# # 
# # clusters_heps <- as.data.frame(kpod[["cluster"]])%>%
# #   rownames_to_column("cell_ID")
# # colnames(clusters_heps)[2] = "cluster"
# 
# ## Export clusters for ML training
# 
# classes = 5
# 
# data.frame(cell_ID = p_heps$yvars, pc1 = p_heps$rotated$PC1, pc2 = p_heps$rotated$PC2) %>%
#   mutate(range = cut_interval(pc1, n = classes))  -> p_bins_5_tmp
# 
# p_bins_5_tmp %>%
#   distinct(range) %>%
#   arrange(range) %>%
#   mutate(bin = c(1:classes)) %>%
#   right_join(p_bins_5) %>%
#   #column_to_rownames("cell_ID") %>%
#   mutate(bin = abs(bin - (classes + 1))) %>%
#   dplyr::rename(cluster = bin) -> p_bins_5
# 
# 
# test <- read_csv("../data/old_stuff/_221115_scDVP_5_clusters.csv")

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
write_csv(meta_heps_cluster, file = "../output/Tables/kmeans_5clusters.csv")
