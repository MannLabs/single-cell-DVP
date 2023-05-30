###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 3E -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/limma_bins_allproteins.R")

## Define number of classes
classes =20

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(-ratio, n = classes)) -> meta_distances_bin

SA_incl_heps <- unique(d$cell_ID)

meta_distances_bin %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bin) %>%
  mutate(bin = abs(bin - (classes + 1))) -> meta_distances_bin

limma_bins_allproteins %>%
  mutate(direction = logFC > 0) %>%
  group_by(direction) %>%
  slice_min(adj.P.Val, n = 10) %>%
  arrange(adj.P.Val) %>%
  drop_na(Symbol) %>%
  dplyr::select(Symbol, adj.P.Val) %>%
  mutate(rank = c(n():1)) -> proteome_top_10

d %>%
  dplyr::select(cell_ID, Symbol, int_core) %>%
  filter(Symbol %in% proteome_top_10$Symbol) %>%
  spread(Symbol, int_core) -> d_top_10

# Plotting expression ratios against true distance
meta_distances_bin %>%
  left_join(d_top_10) %>%
  gather(Symbol, int, (ncol(.)-nrow(proteome_top_10)+1):ncol(.)) %>%
  left_join(limma_bins_allproteins %>% dplyr::select(Symbol, logFC)) %>%
  mutate(zonated = ifelse(logFC > 0, "PV", "CV")) %>%
  replace(is.na(.), 0) %>%
  mutate(int = 2^int) %>%
  group_by(Symbol, bin, zonated) %>%
  summarise(median = median(int), sd = sd(int, na.rm = T)) %>%
  group_by(Symbol) %>%
  left_join(proteome_top_10) %>%
  mutate(sum = sum(median)) %>%
  mutate(ratio = median/sum) -> meta_distances_bin_summary
# group_by(bin, zonated) %>%
#   summarise(ratio_gp = median(ratio), sd_gp = sd(ratio, na.rm = T)) -> meta_distances_bin_summary

meta_distances_bin_summary %>%
  ggplot(aes(x = as.factor(bin), y = ratio, group = Symbol, color = zonated, alpha = rank))+
  geom_point(size = 2)+
  geom_line()+
  #geom_errorbar(aes(ymin = ratio_gp - sd_gp, ymax = ratio_gp + sd_gp), width=.2)+
  scale_color_manual(values = viridis(4)[2:3]) +
  theme_bw()+
  #scale_y_continuous(limits = c(0,0.15)) +
  theme_classic() +
  scale_alpha_continuous(range = c(0,1)) -> plot_expression_top10

ggsave(plot_expression_top10, file = "../Output/Figures/Spatial_expression_top10.pdf", width = 6, height = 5)