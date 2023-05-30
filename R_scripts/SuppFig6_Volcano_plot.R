###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure S6 -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

## Binning

classes = 20

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(-ratio, n = classes)) -> meta_distances_bin

SA_incl_heps <- unique(d$cell_ID)

meta_distances_bin %>%
  filter(cell_ID %in% SA_incl_all) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bin) %>%
  mutate(bin = abs(bin - (classes + 1))) -> meta_distances_bin
  #filter(bin %in% c(1,4)) -> meta_distances_bin

## Limma

d %>%
  filter(cell_ID %in% meta_distances_bin$cell_ID) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(meta_distances_bin$cell_ID)) %>%
  filter(completeness >= 0.5) %>%
  pull(Protein) -> proteome_50_heps

d %>%
  filter(Protein %in% proteome_50_heps) %>%
  filter(cell_ID %in% meta_distances_bin$cell_ID) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide_50

design <- model.matrix(~(meta_distances_bin %>% column_to_rownames("cell_ID"))[colnames(d_wide_50),]$bin)

fit <- lmFit(d_wide_50, design)
fit <- eBayes(fit)
limma <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg) %>%
  arrange(logFC) %>%
  mutate(FC_rank = c(1:nrow(.))) %>%
  mutate(significant = adj.P.Val < 0.05)

## Volcano plotting

ggplot(data = limma, aes(x = logFC, y = -log10(adj.P.Val), fill = -log10(adj.P.Val)))+
  geom_point(alpha = 0.8, pch = 21)+
  theme_classic()+
  scale_fill_viridis(option = "inferno")+
  #geom_hline(yintercept = -log10(0.05), lty = "dotted")+
  geom_vline(xintercept = 0, lty = "dotted") +
  geom_text_repel(data = limma %>% slice_min(adj.P.Val, n = 40),
                  aes(x = logFC, y = -log10(adj.P.Val), label = Symbol), color = "black") -> plot_volcano

ggsave(plot_volcano, file = "../output/Figures/Volcano.pdf", width = 7, height = 6)
