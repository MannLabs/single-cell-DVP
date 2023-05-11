###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances_bins.R")
load("../output/Variables/meta_pg.R")

## Define number of classes
classes = 8

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

d %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_heps)) %>%
  pull(Protein) -> proteome_heps

d %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  filter(Protein %in% proteome_heps) %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide

## Limma statistics
design <- model.matrix(~meta_distances_bins[colnames(d_wide),]$bin)

fit <- lmFit(d_wide, design)
fit <- eBayes(fit)
limma_8bins_allproteins <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

## Save variables
save(limma_8bins_allproteins, file = "../output/Variables/limma_8bins_allproteins.R")
