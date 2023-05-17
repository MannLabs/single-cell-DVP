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
source("../R_scripts/Functions/scale_dplyr.R")
load("../output/Variables/SA_incl_all.R")

img_meta <- read_csv("../data/meta/meta_img-proteome.csv") %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))


## Define number of classes
classes = 20

## Subset to 90% complete proteins
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
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> meta_distances_bins

#save(meta_distances_bins, file = "../output/Variables/meta_distances_bins.R")

d %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_heps)) %>%
  filter(completeness >= 0.7) %>%
  pull(Protein) -> proteome_90_heps

d %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  filter(Protein %in% proteome_90_heps) %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide_90

## Limma statistics
design <- model.matrix(~meta_distances_bins[colnames(d_wide_90),]$bin)

fit <- lmFit(d_wide_90, design)
fit <- eBayes(fit)
limma_8bins_90complete <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg) %>%
  arrange(logFC) %>%
  mutate(FC_rank = c(1:nrow(.))) %>%
  mutate(significant = adj.P.Val < 0.05)

d %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_distances_bins %>% rownames_to_column("cell_ID") %>% drop_na(bin)) %>%
  group_by(Protein, bin) %>%
  summarise(int = log2(median(2^int_core))) %>%
  drop_na(bin) %>%
  spread(bin, int) %>%
  column_to_rownames("Protein") -> d_wide_90

## Heatmappping data
d_heatmap <- d_wide_90[limma_8bins_90complete %>% arrange(logFC) %>% pull(Protein),]

myBreaks <- c(seq(-2,2, by = 0.1))
myColor <- colorRampPalette(viridis(100, option = "inferno"))(length(myBreaks))

pheatmap(scale(t(d_heatmap)),
         breaks = myBreaks,
         color = myColor,
         #cutree_cols = 1,
         show_colnames = T,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F) -> plot_pheatmap

pdf(file = "../Output/Figures/Heatmap_liver-zonation.pdf", width = 8, height = 12)
plot_pheatmap
dev.off()

## Save variables
save(limma_8bins_90complete, file = "../output/Variables/limma_8bins_90complete.R")
save(proteome_90_heps, file = "../output/Variables/proteome_90_heps.R")

## -- Write tables
write_tsv(meta_distances_bins %>%
            rownames_to_column("cell_ID") %>%
            right_join(meta_distances) %>%
            right_join(img_meta) %>%
            mutate(included = cell_ID %in% SA_incl_all) %>%
            mutate(Hepatocyte = cell_ID %in% SA_incl_heps) %>%
            mutate(bin = abs(bin-9)), file = "../output/Tables/scDVP_meta-distances-bin.tsv")
