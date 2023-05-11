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

read_csv("../data/Mouse.MitoCarta3.0.csv") -> mitocarta
  
## OXPHOS
mitocarta %>%
  filter(grepl("^OXPHOS > Complex |^Mitochondrial central dogma > Translation > Mitochondrial ribosome", MitoCarta3.0_MitoPathways)) %>%
  filter(grepl("subunits|Translation", MitoCarta3.0_MitoPathways)) %>%
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

meta_distances_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(d_OXPHOS) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin) %>%
  summarise(median = median(int)) -> d_OXPHOS_tmp

d_OXPHOS_tmp %>%
  filter(bin == 4 | bin == 5) %>%
  group_by(Symbol) %>%
  summarise(int_ref = median(median)) %>%
  right_join(d_OXPHOS_tmp) %>%
  mutate(int_relative = log2(median / int_ref)) %>%
  dplyr::select(bin, Symbol, int_relative) %>%
  drop_na(int_relative) %>%
  spread(Symbol, int_relative) %>%
  column_to_rownames("bin") -> d_OXPHOS

myBreaks <- c(seq(-1,1, by = 0.01))
myColor <- colorRampPalette(viridis(100, option = "inferno"))(length(myBreaks))

as.data.frame(t(scale(d_OXPHOS))) %>%
  rownames_to_column("Protein") %>%
  gather(Zone, zscore, !Protein) %>%
  mutate(Complex = str_replace_all(Protein, "_.*", "")) %>%
  group_by(Complex, Zone) %>%
  summarise(z_median = median(zscore, na.rm = T)) %>%
  ggplot(aes(x = Zone, y = z_median, group = Complex, color = Complex))+
  geom_hline(yintercept = 0)+
  geom_point(size = 4) +
  geom_line(size = 2)+
  scale_color_manual(values = viridis(8)[2:7])+
  theme_classic()

pheatmap(t(scale(d_OXPHOS)), cluster_cols = F, cluster_rows = F, color = myColor, breaks = myBreaks, cellwidth = 10, cellheight = 10) -> plot_OXPHOS
ggsave(plot_OXPHOS, file = "../output/Figures/OXPHOS_hits.pdf", width = 5, height = 20)

##### Alternative approach #####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

load("../output/variables/d.R")
load("../output/Variables/meta_distances_bins.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

read_csv("../data/Mouse.MitoCarta3.0.csv") -> mitocarta

## OXPHOS
mitocarta %>%
  filter(grepl("^OXPHOS > Complex |^Protein import, sorting and homeostasis", MitoCarta3.0_MitoPathways)) %>%
  filter(grepl("subunits|homeostasis", MitoCarta3.0_MitoPathways)) %>%
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

meta_distances_bins %>%
  rownames_to_column("cell_ID") %>%
  left_join(d_OXPHOS) %>%
  mutate(int = 2^int_core) %>%
  group_by(Symbol, bin) %>%
  summarise(median = median(int)) %>%
  filter(median != 1) -> d_OXPHOS_tmp

d_OXPHOS_tmp %>%
  group_by(Symbol) %>%
  summarise(sum_int = sum(median)) %>%
  ungroup() %>%
  right_join(d_OXPHOS_tmp) %>%
  mutate(proportion = median / sum_int) %>%
  mutate(Complex = str_replace_all(Symbol, "_.*", "")) %>%
  group_by(Complex, bin) %>%
  summarise(median_prop = median (proportion, na.rm = T), sd_prop = sd(proportion, na.rm = T)) -> d_OXPHOS

d_OXPHOS %>%
  ggplot(aes(x = bin, y = median_prop, group = Complex, color = Complex))+
  #geom_errorbar(aes(ymin = median_prop - sd_prop, ymax = median_prop + sd_prop), width=.2)+
  geom_hline(yintercept = 0)+
  geom_point(size = 4) +
  geom_line(size = 2)+
  scale_color_manual(values = viridis(8)[2:7])+
  theme_classic()
  