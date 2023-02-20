###########################
#### scDVP Figure Code ####
###########################

#### -- Data wrangling -- ####

## Libraries
library(tidyverse)

## Read relevant data
load("./output/d_all.R")
load("./output/SA_incl_all.R")
load("./output/meta_pg.R")
load("./output/statTable_all.R")

## Additional function
normalize_core <- function(data, proteome_core, SA_incl){
  
  data %>%
    filter(label != "Dimethyl-n-0") %>%
    filter(cell_ID %in% SA_incl) %>%
    drop_na(int) %>%
    filter(Protein %in% proteome_core) %>%
    group_by(cell_ID) %>%
    summarise(median = median(int)) %>%
    ungroup() %>%
    summarise(median_total = median(median)) %>%
    pull(median_total) -> proteome_core_median
  
  data  %>%
    filter(label != "Dimethyl-n-0") %>%
    filter(cell_ID %in% SA_incl) %>%
    drop_na(int) %>%
    filter(Protein %in% proteome_core) %>%
    group_by(cell_ID) %>%
    summarise(median_sample = median(int)) %>%
    ungroup() %>%
    mutate(median_total = proteome_core_median) %>%
    mutate(norm_factor = median_total / median_sample) %>%
    dplyr::select(cell_ID, norm_factor) -> proteome_core_norms
  
  data %>%
    filter(label != "Dimethyl-n-0") %>%
    filter(cell_ID %in% SA_incl) %>%
    drop_na(int) %>%
    left_join(proteome_core_norms) %>%
    mutate(int_core = int * norm_factor) %>%
    dplyr::select(-int, - int) %>%
    mutate(int_core = log2(int_core)) -> d_norm
  
  return(d_norm)
}

## Normalizing hepatocyte data
d_heps <- d_all %>%
  filter(heps == TRUE)

d_heps %>% 
  filter(label != "Dimethyl-n-0") %>%
  filter(cell_ID %in% SA_incl_all) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_heps

d <- normalize_core(data = d_heps, proteome_core = proteome_complete_heps, SA_incl = SA_incl_all) %>%
  left_join(meta_pg)

SA_incl_heps <- unique(d$cell_ID)

d %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("Protein") -> d_complete_heps

meta_heps <- statTable_all %>%
  distinct(bio_ID, label, cell_ID, run_ID, heps) %>%
  filter(cell_ID %in% colnames(d_complete_heps)) %>%
  arrange(cell_ID) %>%
  mutate(run_number = as.numeric(str_replace(run_ID, ".*_", ""))) %>%
  column_to_rownames("cell_ID") 

## Plotting functions

# - Loadings on PC1/PC2
p_heps <- PCAtools::pca(d_complete_heps, metadata = meta_heps, removeVar = 0.1)

PCAtools::biplot(p_heps ,
                 colby = 'heps',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = T,
                 shape = 'heps') +
  scale_color_manual(values = c("grey50","grey50")) -> plot_pca_loadings

ggsave(plot_pca_loadings, file = "./Output/Figures/PCA_Loadings.pdf", width = 7, height = 6)

# - Expression profiles in PCAs
markers <- data.frame(Symbol = c("Ass1", "Cyp2e1"),
                      type = c("PV", "CV")) %>%
           left_join(meta_pg)

d %>%
  left_join(markers) %>%
  drop_na(type) %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  left_join(data.frame(cell_ID = p_heps$yvars, pc1 = p_heps$rotated$PC1)) %>%
  spread(Symbol, int_core) %>%
  arrange(pc1) %>%
  dplyr::select(-pc1) %>%
  column_to_rownames("cell_ID") -> d_pca

meta_heps %>%
  rownames_to_column("cell_ID") %>%
  left_join(d_pca %>% rownames_to_column("cell_ID")) %>%
  column_to_rownames("cell_ID") -> meta_heps_markers

p_heps <- PCAtools::pca(d_complete_heps, metadata = meta_heps_markers, removeVar = 0.1)

# -- Expression profile of Ass1
plot_pca_Ass1 <- PCAtools::biplot(p_heps ,
                   colby = 'Ass1',
                   hline = 0, vline = 0,
                   labSize = 3,
                   lab = NA,
                   encircle = F,
                   encircleFill = F,
                   showLoadings = F) +
                 scale_color_viridis() 

# -- Expression profile of Cyp2e1
plot_pca_Cyp2e1 <- PCAtools::biplot(p_heps ,
                     colby = 'Cyp2e1',
                     hline = 0, vline = 0,
                     labSize = 3,
                     lab = NA,
                     encircle = F,
                     encircleFill = F,
                     showLoadings = F) +
                   scale_color_viridis_c(limits = c(11, 15), oob = scales::squish)

ggsave(file = "./Output/Figures/PCA-Ass1.pdf", plot_pca_Ass1, width = 7, height = 6)
ggsave(file = "./Output/Figures/PCA-Cyp2e1.pdf", plot_pca_Cyp2e1, width = 7, height = 6)

## save variables to file
save(p_heps, file = "./output/p_heps.R")
save(d_complete_heps, file = "./output/d_complete_heps.R")
save(d, file = "./output/d.R")