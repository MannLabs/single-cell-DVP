###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 3A -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../Output/Variables/p_heps.R")
load("../Output/Variables/meta_distances.R")
load("../Output/Variables/statTable_all.R")
load("../Output/Variables/d_complete_heps.R")

## Plotting functions
classes = 8

data.frame(cell_ID = p_heps$yvars, pc1 = p_heps$rotated$PC1, pc2 = p_heps$rotated$PC2) %>%
  mutate(range = cut_interval(pc1, n = classes))  -> p_bins

p_bins %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(p_bins) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> p_bins

# - Scatterplot PC1 and PC2 versus distance
left_join(p_bins %>% rownames_to_column("cell_ID"), meta_distances) %>%
  dplyr::select(cell_ID, ratio, pc1, pc2) %>%
  gather(pc, eigenvalue, 3:4) %>%
  ggplot(aes(x = ratio, y = eigenvalue, color = pc)) +
  #geom_density_2d_filled(contour_var = "count")+
  geom_point(alpha = 0.8, size = 2)+
  geom_smooth(se = FALSE, size = 1, color = "black")+
  theme_classic()+
  labs(x = "Relative distance to PV", y = "PC1 Eigenvalue") +
  facet_wrap(.~pc) +
  scale_color_manual(values = viridis(4)[2:3]) -> plot_distance_pc

ggsave(plot_distance_pc, file = "../Output/Figures/PC1-2_Scatter_GeometricDistance.pdf", width = 6, height = 3)

# - Distance profile in PCA
meta_heps <- statTable_all %>%
  distinct(bio_ID, label, cell_ID, run_ID, heps) %>%
  filter(cell_ID %in% colnames(d_complete_heps)) %>%
  arrange(cell_ID) %>%
  left_join(meta_distances) %>%
  mutate(run_number = as.numeric(str_replace(run_ID, ".*_", ""))) %>%
  drop_na(ratio) %>%
  column_to_rownames("cell_ID") 

p <- PCAtools::pca(d_complete_heps[,rownames(meta_heps)], metadata = meta_heps, removeVar = 0.1)

PCAtools::biplot(p,
                 colby = 'ratio',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = F)+
  scale_color_viridis(option = "viridis") -> plot_pca_geometric_distance

ggsave(plot_pca_geometric_distance, file = "../Output/Figures/PCA_GeometricDistance.pdf", width = 7, height = 6)

## save variables to file
save(p_bins, file = "../output/variables/p_bins.R")

## -- Write tables
write_tsv(as.data.frame(p$rotated) %>% rownames_to_column("sample"), file = "../output/Tables/scDVP-PC_samples.tsv")
write_tsv(as.data.frame(p$loadings) %>% rownames_to_column("Protein"), file = "../output/Tables/scDVP-PC_Proteins.tsv")
