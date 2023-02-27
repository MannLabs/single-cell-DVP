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

## Define number of classes
classes = 8

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> p_bins

p_bins %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(p_bins) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> p_bins

d %>%
  filter(cell_ID %in% rownames(p_bins)) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_heps)) %>%
  pull(Protein) -> proteome_heps

d %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  filter(Protein %in% proteome_heps) %>%
  filter(cell_ID %in% rownames(p_bins)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide

## Limma statistics
design <- model.matrix(~p_bins[colnames(d_wide),]$bin)

fit <- lmFit(d_wide, design)
fit <- eBayes(fit)
limma_8bins_allproteins <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

## Save variables
save(limma_8bins_allproteins, file = "../output/Variables/limma_8bins_allproteins.R")
