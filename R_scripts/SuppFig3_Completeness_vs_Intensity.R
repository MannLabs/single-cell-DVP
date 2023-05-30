###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure S3 -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/Variables/statTable_all.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/CIn_lower_all.R")

statTable_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  drop_na(int) %>%
  filter(label != "DimethNter0") %>%
  mutate(int = log(int, 10)) %>%
  group_by(Protein) %>%
  mutate(completeness = 100*n()/length(SA_incl_all)) %>%
  ungroup() %>%
  group_by(Protein, completeness) %>%
  summarise(int_median = median(int)) %>%
  ungroup() %>%
  drop_na(completeness) %>%
  mutate(bin = cut(completeness, seq(min(completeness)-25, max(completeness) + 25, 25))) %>%
  ggplot(aes(x = bin , y = int_median, fill = bin)) +
  geom_boxplot()+
  theme_bw() +
  scale_fill_manual(values = rev(viridis(5))) -> plot_int_complete

ggsave(plot_int_complete, file = "../output/Figures/Completeness_vs_Intensity.pdf", width = 3, height = 5)
