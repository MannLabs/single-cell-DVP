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

hpa <- read_csv("../data/external/aal3321_thul_sm_table_s7.csv") %>%
  filter(`Cell line` == "U-2 OS") %>%
  dplyr::rename(`Ensembl_human` = `Ensembl id`, Localisation = `IF main protein location`) %>%
  dplyr::select(`Ensembl_human`, Localisation)

orthologs <- read_csv("../data/external/diopt_ortholog-mapping.csv") %>%
  dplyr::rename(Protein = `Search Term`, `Ensembl_human` = `Ensmbl ID  (link HPA)`) %>%
  dplyr::select(Protein, `Ensembl_human`) %>%
  left_join(hpa) %>%
  drop_na(Localisation) %>%
  distinct(Protein, Localisation)

hpa %>%
  mutate(Localisation = ifelse(grepl(";", Localisation), "Other", Localisation)) %>%
  group_by(Localisation) %>%
  summarise(n = n()) %>%
  slice_max(n, n = 10) %>%
  pull(Localisation) -> top_10_hpa

hpa %>%
  mutate(Localisation = ifelse(grepl(";", Localisation), "Other", Localisation)) %>%
  mutate(COI = Localisation %in% top_10_hpa) %>%
  mutate(Localisation = ifelse(COI == FALSE, "Other", Localisation)) %>%
  group_by(Localisation) %>%
  summarise(n = n(), ratio = n / nrow(.)) %>%
  mutate(cat = "HPA") -> ratio_hpa

orthologs %>%
  mutate(Localisation = ifelse(grepl(";", Localisation), "Other", Localisation)) %>%
  mutate(COI = Localisation %in% top_10_hpa) %>%
  mutate(Localisation = ifelse(COI == FALSE, "Other", Localisation)) %>%
  group_by(Localisation) %>%
  summarise(n = n(), ratio = n / nrow(.)) %>%
  mutate(cat = "scDVP") -> ratio_scDVP

ggplot(rbind(ratio_scDVP, ratio_hpa), aes(x = cat, y = ratio, fill = fct_reorder(Localisation, ratio, .desc = TRUE))) +
  geom_bar(stat = "identity") +
  labs(x = "Category", y = "Value", fill = "Group") +
  theme_minimal()+
  scale_fill_manual(values = rev(viridis(10, option = "inferno"))) -> plot_subcellular_bar

ggsave(plot_subcellular_bar, file = "../output/Figures/Subcellular_HPA_scDVP.pdf", width = 5, height = 5)

## Subcellular localisation of significant hits

load("../output/Variables/limma_8bins_allproteins.R")

limma_8bins_allproteins %>%
  left_join(orthologs) %>%
  drop_na(Localisation) %>%
  filter(Localisation %in% top_10_hpa) %>%
  ggplot(aes(x = Localisation, y = B, fill = Localisation))+
  geom_hline(yintercept = 0, lty = "dotted") + 
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = viridis(12)[2:11]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) -> plot_effect_size_by_localisation

ggsave(plot_effect_size_by_localisation, file = "../output/Figures/Subcellular_effect_size.pdf", width = 7, height = 5)
  
  


