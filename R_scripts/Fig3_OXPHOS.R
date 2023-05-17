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

read_csv("../data/external/Mouse.MitoCarta3.0.csv") -> mitocarta

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
  
## OXPHOS
mitocarta %>%
  filter(grepl("^OXPHOS > Complex|^Metabolism > Lipid metabolism", MitoCarta3.0_MitoPathways)) %>%
  filter(grepl("subunits|Lipid", MitoCarta3.0_MitoPathways)) %>%
  mutate(Complex = str_replace_all(str_replace_all(MitoCarta3.0_MitoPathways, "^OXPHOS > ", ""), " > .*", "")) %>%
  dplyr::rename(ENSEMBL = EnsemblGeneID) %>%
  dplyr::select(ENSEMBL, Complex) -> mitocarta_oxphos
  
d %>%
  dplyr::select(cell_ID, ENSEMBL, Symbol, int_core) %>%
  left_join(mitocarta_oxphos) %>%
  drop_na(Complex) %>%
  mutate(Symbol = paste(Complex, Symbol, sep = "_")) %>%
  dplyr::select(-Complex, -ENSEMBL) %>%
  spread(cell_ID, int_core) %>%
  gather(cell_ID, int_core, !Symbol) %>%
  mutate_all(~replace(., is.na(.), 0)) -> d_OXPHOS

meta_distances_bin %>%
  left_join(d_OXPHOS)%>%
  mutate(complex = str_replace_all(Symbol, "_.*", "")) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin, complex) %>%
  summarise(median = median(int), sd = sd(int, na.rm = T)) %>%
  group_by(Symbol) %>%
  mutate(sum = sum(median)) %>%
  mutate(ratio = median/sum) %>%
  filter(sum != classes) %>%
  group_by(bin, complex) %>%
  summarise(ratio_gp = median(ratio), sd_gp = sd(ratio, na.rm = T)) -> meta_distances_bin_summary

meta_distances_bin_summary %>%
  filter(complex == "Complex I") %>%
  mutate(base_CI = ratio_gp) %>%
  dplyr::select(bin, base_CI) %>%
  right_join(meta_distances_bin_summary) %>%
  mutate(ratio_to_CI = ratio_gp / base_CI) -> summary_to_CI


meta_distances_bin_summary %>%
  ggplot(aes(x = as.factor(bin), y = ratio_gp, group = complex, color = complex))+
  geom_point(size = 2)+
  geom_line()+
  #geom_errorbar(aes(ymin = ratio_gp - sd_gp, ymax = ratio_gp + sd_gp), width=.2)+
  scale_color_manual(values = viridis(6)) +
  theme_bw()+
  #scale_y_continuous(limits = c(0,0.15)) +
  theme_classic() +
  scale_alpha_continuous(range = c(0,1)) +
  geom_hline(yintercept = 1/20) -> plot_expression_OXPHOS

## -- Save plots
ggsave(plot_expression_OXPHOS, file = "../output/Figures/OXPHOS_spatial.pdf", width = 6, height = 5)

summary_to_CI %>%
  ggplot(aes(x = as.factor(bin), y = log10(ratio_to_CI), group = complex, color = complex))+
  geom_point(size = 2)+
  geom_line()+
  #geom_errorbar(aes(ymin = ratio_gp - sd_gp, ymax = ratio_gp + sd_gp), width=.2)+
  scale_color_manual(values = viridis(5)) +
  theme_bw()+
  #scale_y_continuous(limits = c(0,0.15)) +
  theme_classic() +
  scale_alpha_continuous(range = c(0,1)) #-> plot_expression_top10


## -- How many mitochondrial proteins covered
d %>%
  drop_na(ENSEMBL) %>%
  distinct(ENSEMBL) %>%
  pull(ENSEMBL) -> all_prots

table(all_prots %in% mitocarta$EnsemblGeneID)
