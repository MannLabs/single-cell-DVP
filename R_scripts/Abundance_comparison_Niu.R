###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/Variables/statTable_all.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

statTable_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  drop_na(int) %>%
  filter(label != "DimethNter0") %>%
  mutate(int = log(int, 10)) %>%
  group_by(Protein) %>%
  summarise(`LFQ intensity Mouse` = mean(int))  -> d_rank_scDVP

d_niu <- read_csv("../data/Niu_msb202210947-sup-0003-datasetev1.csv") %>%
  dplyr::select(`Protein ID`, hHEP) %>%
  mutate(`Ensembl_human` = mapIds(org.Hs.eg.db,
                          keys=str_replace_all(`Protein ID`, ".*;", ""),
                          column="ENSEMBL",
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  mutate(`LFQ intensity hHEP` = log10(hHEP)) %>%
  dplyr::select(-hHEP, -`Protein ID`) %>%
  drop_na(`Ensembl_human`)

read_csv("../data/diopt_ortholog-mapping.csv") %>%
  dplyr::rename(Protein = `Search Term`, `Ensembl_human` = `Ensmbl ID  (link HPA)`) %>%
  dplyr::select(Protein, `Ensembl_human`) %>%
  full_join(d_rank_scDVP) %>%
  drop_na(`Ensembl_human`) %>%
  drop_na(`LFQ intensity Mouse`)%>%
  full_join(d_niu) %>%
  mutate_all(~replace(., is.na(.), 0)) -> d_plot_scatter

cor(d_plot_scatter$`LFQ intensity Mouse`, d_plot_scatter$`LFQ intensity hHEP`, use = "pairwise.complete.obs")

  ggplot() +
    geom_point(data = d_plot_scatter, aes(x = `LFQ intensity Mouse`, y = `LFQ intensity hHEP`)) +
    geom_point(data = (d_plot_scatter %>% mutate_all(~replace(., . == 0, NA))), aes(x = median(`LFQ intensity Mouse`, na.rm = T), y = median(`LFQ intensity hHEP`, na.rm = T)), size = 5, color = "red") +
    theme_classic()+
    geom_abline(slope = 1, intercept = 0, lty = "dotted")+
    scale_y_continuous(breaks = seq(0,12, by = 1)) +
    scale_x_continuous(breaks = seq(0,5, by = 1))+
    geom_hline(yintercept = median(d_plot_scatter %>% mutate_all(~replace(., . == 0, NA)) %>% pull(`LFQ intensity hHEP`), na.rm = T)) +
    geom_vline(xintercept = median(d_plot_scatter %>% mutate_all(~replace(., . == 0, NA)) %>% pull(`LFQ intensity Mouse`), na.rm = T)) -> plot_scatter
  
read_csv("../data/diopt_ortholog-mapping.csv") %>%
    dplyr::rename(Protein = `Search Term`, `Ensembl_human` = `Ensmbl ID  (link HPA)`) %>%
    dplyr::select(Protein, `Ensembl_human`) %>%
    full_join(d_rank_scDVP) %>%
    drop_na(`Ensembl_human`) %>%
    drop_na(`LFQ intensity Mouse`)%>%
    full_join(d_niu) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    gather(dataset, int, 3:4) %>%
    filter(int != 0) %>%
    ggplot(aes(x = dataset, y = int)) +
    geom_boxplot()+
    theme_classic() +
    scale_y_continuous(breaks = seq(0,12, by = 1))

read_csv("../data/diopt_ortholog-mapping.csv") %>%
  dplyr::rename(Protein = `Search Term`, `Ensembl_human` = `Ensmbl ID  (link HPA)`) %>%
  dplyr::select(Protein, `Ensembl_human`) %>%
  full_join(d_rank_scDVP) %>%
  drop_na(`Ensembl_human`) %>%
  drop_na(`LFQ intensity Mouse`)%>%
  full_join(d_niu) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(Human_Mouseis0 = `LFQ intensity Mouse` == 0) %>%
  filter(`LFQ intensity hHEP` != 0) %>%
  ggplot(aes(x = Human_Mouseis0, y = `LFQ intensity hHEP`, fill = Human_Mouseis0)) +
  geom_boxplot()+
  theme_classic() +
  scale_y_continuous(breaks = seq(0,12, by = 1)) +
  scale_fill_manual(values = viridis(4)[2:3]) +
  labs(title = "Intensity of hHep proteins when scDVP protein not detected") -> plot_Human_Mouseis0

read_csv("../data/diopt_ortholog-mapping.csv") %>%
  dplyr::rename(Protein = `Search Term`, `Ensembl_human` = `Ensmbl ID  (link HPA)`) %>%
  dplyr::select(Protein, `Ensembl_human`) %>%
  full_join(d_rank_scDVP) %>%
  drop_na(`Ensembl_human`) %>%
  drop_na(`LFQ intensity Mouse`)%>%
  full_join(d_niu) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(Mouse_Humanis0 = `LFQ intensity hHEP` == 0) %>%
  filter(`LFQ intensity Mouse` != 0) %>%
  ggplot(aes(x = Mouse_Humanis0, y = `LFQ intensity Mouse`, fill = Mouse_Humanis0)) +
  geom_boxplot()+
  theme_classic() +
  scale_y_continuous(breaks = seq(0,12, by = 1)) +
  scale_fill_manual(values = viridis(4)[2:3]) +
  labs(title = "Intensity of scDVP proteins when hHep protein not detected") -> plot_Mouse_Humanis0
  
ggsave(plot_scatter, file = "../output/Figures/Abundance_comparison_scatter.pdf", height = 5, width = 5)
ggsave(plot_Human_Mouseis0, file = "../output/Figures/Abundance_comparison_HumanUnique.pdf", height = 5, width = 3)
ggsave(plot_Mouse_Humanis0, file = "../output/Figures/Abundance_comparison_scDVPunique.pdf", height = 5, width = 3)
