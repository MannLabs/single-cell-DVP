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
load("../output/Variables/SA_incl_all.R")

## Define number of classes
classes = 20

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
  filter(completeness >= 0.9) %>%
  pull(Protein) -> proteome_90_heps

d %>%
  filter(cell_ID %in% SA_incl_all) %>%
  filter(Protein %in% proteome_90_heps) %>%
  filter(cell_ID %in% rownames(p_bins)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide_90

## Limma statistics
design <- model.matrix(~p_bins[colnames(d_wide_90),]$bin)

fit <- lmFit(d_wide_90, design)
fit <- eBayes(fit)
limma_8 <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg) %>%
  arrange(logFC) %>%
  mutate(FC_rank = c(1:nrow(.))) %>%
  mutate(significant = adj.P.Val < 0.05)

d %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(p_bins %>% rownames_to_column("cell_ID") %>% drop_na(bin)) %>%
  group_by(Protein, bin) %>%
  summarise(int = log2(median(2^int_core))) %>%
  drop_na(bin) %>%
  spread(bin, int) %>%
  column_to_rownames("Protein") -> d_wide_90

## Heatmappping data
d_heatmap <- d_wide_90[limma_8 %>% arrange(logFC) %>% pull(Protein),]

myBreaks <- c(seq(-2,2, by = 0.1))
myColor <- colorRampPalette(viridis(100, option = "inferno"))(length(myBreaks))

pheatmap(scale(t(d_heatmap)),
         breaks = myBreaks,
         color = myColor,
         cutree_cols = 1,
         show_colnames = T,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F) -> plot_pheatmap

ggsave(plot_pheatmap, file = "../Output/Figures/Heatmap_liver-zonation.pdf", width = 8, height = 12)