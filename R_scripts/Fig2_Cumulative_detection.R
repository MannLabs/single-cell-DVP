###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 2A -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/statTable_all.R")
load("../output/variables/SA_incl_all.R")
load("../output/variables/d.R")

statTable_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  group_by(Protein) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  mutate(rank = c(nrow(.):1)) %>%
  mutate(n = n/max(n)) %>%
  ggplot(aes(x = rank, y = n))+
  geom_point()+
  theme_classic() -> plot

ggsave(plot, file = "../output/Figures/Rank_versus_completeness.pdf", width = 5, height = 5)

statTable_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  group_by(Protein) %>%
  summarise(n = n(), int_median =log2(median(int))) %>%
  ggplot(aes(x = int_median, y = n))+
  geom_point()

## Contamination with blood

#HBA P01942
#HBB P01942
#FGA E9PV24 !
Fgb
Fgg

d %>%
  filter(cell_ID %in% SA_incl_all) %>%
  filter(Symbol %in% c("Hba", "Hbb", "Hbc", "Hbd", "Fga", "Fgb", "Fgg")) %>%
  ggplot(aes(x = cell_ID, y = int_core, color = Symbol))+
  geom_point()

statTable_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  drop_na(int) %>%
  filter(label != "DimethNter0") %>%
  mutate(POI = Symbol %in% c("Hba", "Hbb", "Hbc", "Hbd", "Fga", "Fgb", "Fgg", "Alb")) %>%
  mutate(int = log(int, 10)) %>%
  group_by(Protein, POI, Symbol) %>%
  summarise(`log10 LFQ intensity` = mean(int)) %>%
  ungroup() %>%
  arrange(`log10 LFQ intensity`) %>%
  mutate(`Protein rank` = c(length(unique(Protein)):1)) -> fig_ProteinRank

ggplot()+
  geom_point(data = fig_ProteinRank %>% filter(POI == FALSE), aes(x = `Protein rank`, y = `log10 LFQ intensity`), color = "grey20", alpha = 0.6)+
  geom_point(data = fig_ProteinRank %>% filter(POI == TRUE), aes(x = `Protein rank`, y = `log10 LFQ intensity`), color = "orange", alpha = 1, size = 3) +
  geom_label_repel(data = (fig_ProteinRank %>% filter(POI == TRUE)), aes(x = `Protein rank`, y = `log10 LFQ intensity`, label = Symbol),max.overlaps = 20) +
  theme_bw() -> plot_blood

ggsave(plot_blood, file = "../output/Figures/Blood_contamination.pdf", width = 3, height = 3)






