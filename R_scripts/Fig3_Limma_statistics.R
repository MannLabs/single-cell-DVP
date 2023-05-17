###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/meta_pg.R")

classes = 20

SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> meta_distances_bins

meta_distances_bins %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bins) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") -> meta_distances_bins

## Filter on 70% complete proteins 
d %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  column_to_rownames("Protein") %>%
  filter(rowSums(is.na(.)) / ncol(.) <= 0.5) -> d_wide_70

## Limma statistics
design <- model.matrix(~meta_distances_bins[colnames(d_wide_70),]$bin)

fit <- lmFit(d_wide_70, design)
fit <- eBayes(fit)
limma_bins_allproteins <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

## Save variables
save(limma_bins_allproteins, file = "../output/Variables/limma_8bins_allproteins.R")

## -- Write tables
write_tsv(limma_bins_allproteins, file = "../output/Tables/scDVP_ANOVA.tsv")
