###########################
#### scDVP Figure Code ####
###########################

## -- Figure 3G, Supplementary Figure S8 -- ##

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/SA_incl_all.R")
load("../output/Variables/meta_pg.R")

## Data binning
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

mitocarta_loc <- read_csv("../data/external/Mouse.MitoCarta3.0_subcellular.csv") %>%
  dplyr::select(EnsemblGeneID, `HPA_Main_Location_2020 (Reliability)`) %>%
  dplyr::rename(Compartment = `HPA_Main_Location_2020 (Reliability)`, ENSEMBL = EnsemblGeneID) %>%
  filter(grepl("Supported|Approved|Enhanced", Compartment)) %>%
  mutate(Compartment = str_replace_all(Compartment, " \\(.*", "")) %>%
  separate_rows(Compartment, sep = ";")

mitocarta_loc %>%
  group_by(Compartment) %>%
  summarise(n = n()) %>%
  slice_max(n, n = 10) %>%
  pull(Compartment) -> top_10_compartments 

d %>%
  left_join(mitocarta_loc) %>%
  filter(Compartment %in% top_10_compartments) %>%
  right_join(meta_distances_bin) %>%
  drop_na(Compartment) %>%
  group_by(bin, Compartment) %>%
  summarise(int_sum = sum(int_core)) %>%
  ungroup() %>%
  group_by(bin) %>%
  mutate(int_bin = sum(int_sum)) %>%
  mutate(ratio = int_sum / int_bin) -> annotations_scDVP

annotations_scDVP %>%
  group_by(Compartment) %>%
  summarise(mean = mean(ratio))

library_proteome <- read_tsv("../data/external/scDVP_report.pg_matrix.tsv") %>%
  dplyr::select(-c(2:5)) %>%
  gather(sample, int, !Protein.Group) %>%
  group_by(Protein.Group) %>% 
  summarise(int = median(int, na.rm = T)) %>%
  mutate(ENSEMBL = mapIds(org.Mm.eg.db,
                          keys=str_replace_all(Protein.Group, ".*;", ""),
                          column="ENSEMBL",
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  left_join(mitocarta_loc)

library_proteome %>%
  filter(Compartment %in% top_10_compartments) %>%
  drop_na(Compartment) %>%
  group_by(Compartment) %>%
  summarise(int = sum(int)) %>%
  ungroup() %>%
  mutate(ratio = int / sum(int)) %>%
  mutate(bin = 0) -> annotations_library

annotations_scDVP %>%
  full_join(annotations_library) -> annotations_full

ggplot(annotations_full, aes(x = bin, y = ratio, alpha = fct_reorder(Compartment, ratio, .desc = TRUE))) +
  geom_bar(stat = "identity", fill = viridis(5)[2], color = "grey20") +
  labs(x = "Category", y = "Value", fill = "Group") +
  theme_classic()+
  scale_alpha_manual(values = seq(0,1,by = 0.1)) -> plot_scDVP

ggplot(data = annotations_scDVP %>% group_by(Compartment) %>% mutate(z = scale(ratio)), aes(x = bin, y = z)) + 
  geom_hline(yintercept = 0, lty = "dotted") +
  geom_line()+
  geom_smooth(method = "lm")+
  facet_wrap(.~Compartment, ncol = 5)+
  theme_classic() -> plot_lm

ggplot(data = annotations_scDVP, aes(x = bin, y = ratio)) + 
  geom_line()+
  geom_smooth(method = "lm")+
  facet_wrap(.~Compartment, ncol = 5, scales = "free_y")+
  theme_classic()

annotations_scDVP %>%
  group_by(Compartment) %>%
  mutate(z = scale(ratio)) %>%
  do(model = summary(lm(ratio ~ bin, data = .))) %>%
  mutate(r.squared = model$r.squared) %>% #-> annotations_scDVP_rsqu
  mutate(p.value = model$coefficients[,4][2])

## -- Another alternative approach, map from MouseMitocarta 3.0 ;; applied to library proteome

library_proteome <- read_tsv("../data/external/scDVP_report.pg_matrix.tsv") %>%
  dplyr::select(-c(2:5)) %>%
  gather(sample, int, !Protein.Group) %>%
  group_by(Protein.Group) %>% 
  summarise(int = median(int, na.rm = T)) %>%
  mutate(ENSEMBL = mapIds(org.Mm.eg.db,
                          keys=str_replace_all(Protein.Group, ".*;", ""),
                          column="ENSEMBL",
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  left_join(mitocarta_loc)

library_proteome %>%
  filter(Compartment %in% top_10_compartments) %>%
  drop_na(Compartment) %>%
  group_by(Compartment) %>%
  summarise(int = sum(int)) %>%
  ungroup() %>%
  mutate(ratio = int / sum(int)) -> annotations_library

ggplot(annotations_library, aes(x = 1, y = ratio, fill = fct_reorder(Compartment, ratio, .desc = TRUE))) +
  geom_bar(stat = "identity") +
  labs(x = "Category", y = "Value", fill = "Group") +
  theme_minimal()+
  scale_fill_manual(values = rev(viridis(10, option = "cividis"))) -> plot_library

## Save plots
ggsave(plot_library, file = "../output/Figures/Localisation_library.pdf", width = 3, height = 5)
ggsave(plot_scDVP, file = "../output/Figures/Localisation_scDVP.pdf", width = 11, height = 5)
ggsave(plot_lm, file = "../output/Figures/Localisation_lm.pdf", width = 8, height = 5)

