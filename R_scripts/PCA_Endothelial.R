###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d_all_norm.R")
load("../output/variables/statTable_all.R")

## Endothelial cell PCA
d_all_norm %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("Protein") -> d_complete_all

meta_all <- statTable_all %>%
  distinct(bio_ID, label, cell_ID, run_ID, heps) %>%
  filter(cell_ID %in% colnames(d_complete_all)) %>%
  arrange(cell_ID) %>%
  mutate(run_number = as.numeric(str_replace(run_ID, ".*_", ""))) %>%
  column_to_rownames("cell_ID") 

p <- PCAtools::pca(d_complete_all, metadata = meta_all, removeVar = 0.1)

PCAtools::biplot(p,
                 colby = 'heps',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = F,
                 shape = 'heps', shapekey = c('FALSE'=25, 'TRUE'=21), fill = 'heps') +
          scale_color_manual(values = c("darkred","grey")) -> pca_endothelial

ggsave(pca_endothelial, file = "../Output/Figures/PCA_Endothelial.pdf", width = 7, height = 6)