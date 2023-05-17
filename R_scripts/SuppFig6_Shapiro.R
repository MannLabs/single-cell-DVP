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
load("../output/Variables/proteome_90_heps.R")
load("../output/Variables/limma_8bins_90complete.R")

# Binning
classes = 20

## Subset to 90% complete proteins
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

## Shapiro-Wilk-Test
d %>%
  filter(Protein %in% proteome_90_heps) %>%
  group_by(Protein) %>%
  summarise(shapiro.p = shapiro.test(as.vector(2^int_core))$p.value, shapiro.stat = shapiro.test(as.vector(2^int_core))$statistic) %>%
  mutate(p.adj = p.adjust(shapiro.p, method = "fdr")) %>%
  left_join(limma_8bins_90complete) %>%
  mutate(significant_shapiro = p.adj < 0.05) %>%
  mutate(significant_ANOVA = adj.P.Val < 0.05) %>%
  arrange(desc(shapiro.p)) -> d_shapiro

ggplot(d_shapiro, aes(x = shapiro.stat, y = shapiro.p))+
  geom_point() +
  theme_bw()+
  scale_color_manual(values = viridis(4)[2:3])+
  geom_text_repel(data = d_shapiro[1:10,], aes(x = shapiro.stat, y = shapiro.p, label = Symbol)) -> plot_shapiro

ggsave(plot_shapiro, file = "../output/Figures/Shapiro.pdf", width = 5, height = 5)

d_shapiro %>%
  arrange(-shapiro.stat) %>%
  top_n(shapiro.stat, n = 10) %>%
  pull(Protein) -> shapiro_top10

d %>%
  dplyr::select(cell_ID, Protein, int_core) %>%
  filter(Protein %in% shapiro_top10) %>%
  spread(Protein, int_core) -> d_shapiro_top_10

meta_distances_bin %>%
  left_join(d_shapiro_top_10) %>%
  gather(Protein, int, (ncol(.)-length(shapiro_top10)+1):ncol(.)) %>%
  left_join(limma_8bins_90complete %>% dplyr::select(Protein, logFC)) %>%
  mutate(zonated = ifelse(logFC > 0, "PV", "CV")) %>%
  replace(is.na(.), 0) %>%
  mutate(int = 2^int) %>%
  group_by(Protein, bin, zonated) %>%
  summarise(median = median(int), sd = sd(int)) %>%
  group_by(Protein) %>%
  mutate(sum = sum(median)) %>%
  mutate(ratio = median/sum) %>%
  group_by(bin, zonated) %>%
  summarise(ratio_gp = median(ratio), sd_gp = sd(ratio)) -> meta_distances_bins_shapiro_summary

meta_distances_bins_shapiro_summary %>%
  ggplot(aes(x = as.factor(bin), y = ratio_gp, group = zonated, color = zonated))+
  geom_point(size = 2)+
  geom_line()+
  geom_errorbar(aes(ymin = ratio_gp - sd_gp, ymax = ratio_gp + sd_gp), width=.2)+
  scale_color_manual(values = viridis(4)[2:3]) +
  theme_bw()+
  scale_y_continuous(limits = c(0,0.3)) -> plot_shapiro_top10

ggsave(plot_shapiro_top10, file = "../output/Figures/Shapiro_top10.pdf", width = 7, height = 5)
