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
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

## Fatty acid metabolism
proteome_FA <- c("Acly", "Abcd3", "Acaa1a", "Acaa1b", "Ehhadh", "Hsd17b4", "Hmgcs1", "Hmgcs2", "Hmgcl", "Acat1", "Cat", "Etfb", "Sdhb", "Acox1", "Ephx1")

d %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  filter(Symbol %in% proteome_FA) -> d_FA

meta_distances_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(d_FA) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin) %>%
  summarise(median = median(int)) -> d_peroxisome_tmp

d_peroxisome_tmp %>%
  filter(bin == 4 | bin == 5) %>%
  group_by(Symbol) %>%
  summarise(int_ref = median(median)) %>%
  #dplyr::select(-bin) %>%
  right_join(d_peroxisome_tmp) %>%
  mutate(int_relative = log2(median / int_ref)) %>%
  dplyr::select(bin, Symbol, int_relative) %>%
  spread(Symbol, int_relative) %>%
  column_to_rownames("bin") -> d_peroxisome

myBreaks <- c(seq(-1,1, by = 0.01))
myColor <- colorRampPalette(viridis(100, option = "inferno"))(length(myBreaks))

pheatmap(t(d_peroxisome), cluster_cols = F, cluster_rows = F, color = myColor, breaks = myBreaks, cellwidth = 10, cellheight = 10) -> plot_peroxisomes
ggsave(plot_peroxisomes, file = "../output/Figures/Peroxisomal_hits.pdf")

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances_bins.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

## Urea cycle
proteome_UC <- c("Cps1", "Nags", "Otc", "Ornt1", "Arg1", "Ass1", "Asl", "Fh1", "Mdh", "Mdh2", "Ast", "Oat", "Glul", "Gls2")

d %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  filter(Symbol %in% proteome_UC) -> d_UC

meta_distances_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(d_UC) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin) %>%
  summarise(median = median(int)) -> d_UC_tmp

d_UC_tmp %>%
  filter(bin == 4 | bin == 5) %>%
  group_by(Symbol) %>%
  summarise(int_ref = median(median)) %>%
  #dplyr::select(-bin) %>%
  right_join(d_UC_tmp) %>%
  mutate(int_relative = log2(median / int_ref)) %>%
  dplyr::select(bin, Symbol, int_relative) %>%
  spread(Symbol, int_relative) %>%
  column_to_rownames("bin") -> d_UC

myBreaks <- c(seq(-2,2, by = 0.001))
myColor <- colorRampPalette(viridis(100, option = "inferno"))(length(myBreaks))

pheatmap(t(scale(d_UC)), cluster_cols = F, cluster_rows = F, color = myColor, breaks = myBreaks, cellwidth = 10, cellheight = 10) -> plot_urea
ggsave(plot_urea, file = "../output/Figures/Ureacycle_hits.pdf")
