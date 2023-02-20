###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure XX -- ####
library(tidyverse)

## Read relevant data
load("./output/d_all_norm.R")
load("./output/statTable_all.R")

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

cairo_pdf("./Output/Figures/PCA_Endothelial.pdf", 7, 6)
PCAtools::biplot(p,
                 colby = 'heps',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = F,
                 shape = 'heps', shapekey = c('FALSE'=25, 'TRUE'=21), fill = 'heps') +
          scale_color_manual(values = c("darkred","grey"))
dev.off()