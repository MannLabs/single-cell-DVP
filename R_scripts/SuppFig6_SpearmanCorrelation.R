###########################
#### scDVP Figure Code ####
###########################

#### -- Figure S5 -- ####

## -- Prepare Workspace

cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")

# Binning
classes = 4

SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> meta_distances_bin

meta_distances_bin %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bin) %>%
  filter(cell_ID %in% SA_incl_heps)  -> meta_distances_bin

d %>%
  left_join(meta_distances_bin) %>%
  drop_na(bin) %>%
  group_by(Protein, bin) %>%
  summarise(median = median(int_core)) %>%
  spread(bin, median) %>%
  column_to_rownames("Protein") -> d_wide_median

cormat <- round(cor(d_wide_median, use = "pairwise.complete", method = "spearman"),2)

cairo_pdf("../output/Figures/Zone_correlation_heatmap_viridis.pdf", 5, 5)
pheatmap(cormat, color = viridis(50), show_rownames = T, show_colnames = F, cluster_rows = F, cluster_cols = F,
         cellwidth = 10, cellheight = 10, border_color = "white")
dev.off()
