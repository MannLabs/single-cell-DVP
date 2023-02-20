###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure XX -- ####
library(tidyverse)
library(ggrepel)

## Read relevant data
load("./output/statTable_all.R")
load("./output/SA_incl_all.R")
load("./output/meta_pg.R")

## Plotting functions
statTable_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  drop_na(int) %>%
  filter(label != "DimethNter0") %>%
  mutate(int = log(int, 10)) %>%
  group_by(Protein) %>%
  summarise(`LFQ intensity` = mean(int)) %>%
  ungroup() %>%
  arrange(`LFQ intensity`) %>%
  mutate(`Protein rank` = c(length(unique(Protein)):1)) %>%
  left_join(meta_pg) %>%
  mutate(transcription = grepl("transcription", Genename)) -> fig_ProteinRank

ggplot()+
  geom_point(data = fig_ProteinRank %>% filter(transcription == FALSE), aes(x = `Protein rank`, y = `LFQ intensity`), color = "grey20", alpha = 0.6)+
  geom_point(data = fig_ProteinRank %>% filter(transcription == TRUE), aes(x = `Protein rank`, y = `LFQ intensity`), fill  = viridis(3)[3], size = 4, pch = 21, color= "black")+
  geom_label_repel(data = (fig_ProteinRank %>% arrange(`LFQ intensity`) %>% slice_head(n = 10)), aes(x = `Protein rank`, y = `LFQ intensity`, label = Symbol),max.overlaps = 20) +
  geom_text_repel(data = (fig_ProteinRank %>% filter(transcription == TRUE)), aes(x = `Protein rank`, y = `LFQ intensity`, label = Symbol),max.overlaps = 20) +
  geom_label_repel(data = (fig_ProteinRank %>% arrange(`LFQ intensity`) %>% slice_tail(n = 10)), aes(x = `Protein rank`, y = `LFQ intensity`, label = Symbol),max.overlaps = 20) +
  theme_bw() -> plot_abundanceRange

## Save plot to file
ggsave(plot_abundanceRange, file = "./output/Figures/Rank_versus_Intensity.pdf", width = 5, height = 5)
