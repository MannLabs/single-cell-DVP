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
load("../output/Variables/meta_pg.R")

# Binning
classes = 20

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

## Fatty acid metabolism
proteome_FA <- c("Acly", "Abcd3", "Acaa1a", "Acaa1b", "Ehhadh", "Hsd17b4", "Hmgcs1", "Hmgcs2", "Hmgcl", "Acat1", "Cat", "Etfb", "Sdhb", "Acox1", "Ephx1")

d %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  filter(Symbol %in% proteome_FA) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Symbol) %>%
  replace_na(list(int_core = 0)) -> d_FA

meta_distances_bin %>%
  left_join(d_FA) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin) %>%
  summarise(median = median(int)) -> d_peroxisome_tmp

d_peroxisome_tmp %>%
  filter(bin == 10 | bin == 11) %>%
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

pdf(file = "../output/Figures/Peroxisomal_hits.pdf")
pheatmap(t(d_peroxisome), cluster_cols = F, cluster_rows = F, color = myColor, breaks = myBreaks, cellwidth = 10, cellheight = 10, border_color = "white" ) #-> plot_peroxisomes
dev.off()

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

# Binning
classes = 20

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

## Urea cycle
proteome_UC <- c("Cps1", "Nags", "Otc", "Ornt1", "Arg1", "Ass1", "Asl", "Fh1", "Mdh", "Mdh2", "Ast", "Oat", "Glul", "Gls2")

d %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  filter(Symbol %in% proteome_UC) -> d_UC

d %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  filter(Symbol %in% proteome_UC) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Symbol) %>%
  left_join(meta_distances_bin) %>%
  drop_na(bin) %>%
  group_by(Symbol, bin) %>%
  summarise(proportion_NAs = mean(is.na(int_core))) %>%
  filter(proportion_NAs < 0.5) %>%
  mutate(include_bin = paste(Symbol, bin, sep = "_")) %>%
  pull(include_bin) -> UC_NAs

meta_distances_bin %>%
  left_join(d_UC) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin) %>%
  summarise(median = median(int)) -> d_UC_tmp

d_UC_tmp %>%
  mutate(bin_ID = paste(Symbol, bin, sep = "_")) %>%
  filter(bin_ID %in% UC_NAs) %>%
  dplyr::select(- bin_ID) %>%
  group_by(Symbol) %>%
  summarise(int_ref = median(median, na.rm = T)) %>%
  #dplyr::select(-bin) %>%
  right_join(d_UC_tmp) %>%
  mutate(bin_ID = paste(Symbol, bin, sep = "_")) %>%
  filter(bin_ID %in% UC_NAs) %>%
  dplyr::select(- bin_ID) %>%
  mutate(int_relative = log2(median / int_ref)) %>%
  dplyr::select(bin, Symbol, int_relative) %>%
  spread(Symbol, int_relative) %>%
  column_to_rownames("bin") -> d_UC

myBreaks <- c(seq(-2,2, by = 0.001))
myColor <- colorRampPalette(viridis(100, option = "inferno"))(length(myBreaks))

pdf(file = "../output/Figures/Ureacycle_hits.pdf")
pheatmap(t(scale(d_UC)), cluster_cols = F, cluster_rows = F, color = myColor, breaks = myBreaks, cellwidth = 10, cellheight = 10, border_color = "white") #-> plot_urea
dev.off()




                        