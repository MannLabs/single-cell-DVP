###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d_all.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/meta_pg.R")
load("../output/variables/SA_incl_all.R")
load("../output/variables/img_fluovalues.R")
meta_binuc <- read_csv("../data/meta_binucleation.csv") %>%
  filter(is.na(Notes) & (Classification == "Binuc" | Classification == "Mono"))  %>%
  mutate(cell_ID = str_replace(cell_ID, "Dimethyl-n-", "target"))

## Additional function
source("./Functions/normalize_core_mean.R")

############################
## All with normalization ##
############################

## Define new number of classes

SA_incl_binuc <- d_all %>%
  filter(cell_ID %in% SA_incl_all) %>%
  filter(heps == TRUE) %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  filter(cell_ID %in% meta_binuc$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

classes = 4

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> meta_distances_bins

meta_distances_bins %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bins) %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  mutate(bin = abs(bin - (classes + 1))) %>%
  left_join(meta_binuc) -> meta_distances_bins_4

## Data wrangling and normalization

d_all %>%
  filter(heps == TRUE) %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_binuc

## From Kati's DVP dataset, take the stable core proteome
load(file = "Z:/Florian_Kati/1_Binucleation/1_Analysis/DVP_binuc_proteome-stable.R")

d_all %>%
  group_by(cell_ID) %>%
  summarise(median_int = median(int, na.rm = T)) %>%
  ggplot(aes(x = median_int))+
  geom_histogram()
  
d_all %>%
  group_by(cell_ID) %>%
  summarise(median_int = median(int, na.rm = T)) %>%
  left_join(meta_distances_bins_4) %>%
  drop_na(Classification) %>%
  filter(median_int < 6000 & median_int > 1900) %>%
  pull(cell_ID) -> SA_incl_binuc

d_norm <- normalize_core_mean(data = d_all,
                          proteome_core = proteome_stable,
                          SA_incl = SA_incl_binuc) %>%
  left_join(meta_distances_bins_4)

ggplot(data = d_norm %>% filter(Protein %in% proteome_stable), aes(x = cell_ID, y = int_core, fill = Classification))+
  geom_boxplot() +
  labs(title = "Median intensities after separate normalization, complete proteome")

ggplot(data = d_norm, aes(x = cell_ID, y = int_core, fill = Classification))+
  geom_boxplot() +
  labs(title = "Median intensities after separate normalization, whole proteome")

d_norm %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_binuc)) %>%
  filter(completeness > 0.3) %>%
  pull(Protein) -> proteome_heps_30

d_norm %>%
  left_join(meta_distances_bins_4) %>%
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ Classification,
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_complete

ggplot(data = wilcox_complete, aes(x = estimate, y = -log10(p.adj)) )+
  geom_point()



## Normalization within each group

# Identify shoulder in density
d_all %>%
  group_by(cell_ID) %>%
  summarise(median_int = median(int, na.rm = T)) %>%
  left_join(meta_distances_bins_4) %>%
  drop_na(Classification) %>%
  ggplot(aes(x = median_int, fill = Classification))+
  geom_density(alpha = 0.3)

d_all %>%
  group_by(cell_ID) %>%
  summarise(median_int = median(int, na.rm = T)) %>%
  left_join(meta_distances_bins_4) %>%
  drop_na(Classification) %>%
  filter(Classification == "Mono") %>%
  filter(median_int > 2200 & median_int < 2700) %>%
  pull(cell_ID) -> SA_incl_mono

d_all %>%
  group_by(cell_ID) %>%
  summarise(median_int = median(int, na.rm = T)) %>%
  left_join(meta_distances_bins_4) %>%
  drop_na(Classification) %>%
  filter(Classification == "Binuc") %>%
  filter(median_int < 2900 & median_int > 1900) %>%
  pull(cell_ID) -> SA_incl_binuc

d_all %>%
  filter(cell_ID %in% c(binuc_incl, mono_incl)) %>%
  group_by(cell_ID) %>%
  summarise(median_int = median(int, na.rm = T)) %>%
  left_join(meta_distances_bins_4) %>%
  drop_na(Classification) %>%
  ggplot(aes(x = median_int, fill = Classification))+
  geom_density(alpha = 0.3)

# Normalize in each class
d_all %>%
  filter(heps == TRUE) %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_binuc

d_binuc <- normalize_core_mean(data = d_all,
                          proteome_core = proteome_complete_binuc,
                          SA_incl = SA_incl_binuc) %>%
  left_join(meta_distances_bins_4)

d_all %>%
  filter(heps == TRUE) %>%
  filter(cell_ID %in% SA_incl_mono) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_mono

d_mono <- normalize_core_mean(data = d_all,
                          proteome_core = proteome_complete_binuc,
                          SA_incl = SA_incl_mono) %>%
  left_join(meta_distances_bins_4)

# Join the two datasets

d <- rbind(d_binuc, d_mono) %>%
  dplyr::rename(int = int_core) %>% 
  dplyr::select(-norm_factor)

ggplot(data = d %>% filter(Protein %in% proteome_complete_mono), aes(x = cell_ID, y = int, fill = Classification))+
  geom_boxplot()

d %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_all

d_norm <- normalize_core_mean(data = d %>% mutate(int = 2^int),
                         proteome_core = proteome_complete_all,
                         SA_incl = unique(d$cell_ID))

ggplot(data = d_norm %>% filter(Protein %in% proteome_complete_mono), aes(x = cell_ID, y = int_core, fill = Classification))+
  geom_boxplot() +
  labs(title = "Median intensities after separate normalization, complete proteome")

ggplot(data = d_norm, aes(x = cell_ID, y = int_core, fill = Classification))+
  geom_boxplot() +
  labs(title = "Median intensities after separate normalization, whole proteome")

## PCA
meta_distances_bins_4 %>%
  column_to_rownames("cell_ID") -> meta_pca

d_norm %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  filter(cell_ID %in% rownames(meta_pca)) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("Protein") -> d_complete_norm

# Plotting functions

p_heps <- PCAtools::pca(d_complete_norm, metadata = meta_pca[colnames(d_complete_norm),], removeVar = 0.1)

PCAtools::biplot(p_heps ,
                 colby = 'Classification',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = F,
                 legendPosition = 'right')+
  scale_color_manual(values = viridis(4)[2:3])  -> plot_pca12

d_norm %>%
  filter(Protein %in% proteome_heps_50) %>%
  group_by(Protein) %>%
  t_test(int_core ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_complete

ggplot(wilcox_complete, aes(x = estimate, y = -log10(p.adj)))+
  geom_point()+
  labs(title = "Binuc vs mono of core proteome")

mean(wilcox_complete$estimate)

d_norm %>%
  filter(Protein == "P43276") %>%
  ggplot(aes(x = Classification, y = int_core))+
  geom_boxplot()+
  labs(title = "Levels of marker histone H1f5 (should be elevated in binucleated cells)")
  


















d_binuc <- normalize_core(data = d_all,
                          proteome_core = proteome_complete_binuc,
                          SA_incl = SA_incl_binuc) %>%
  left_join(meta_distances_bins_4)

d_binuc %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_binuc)) %>%
  filter(completeness > 0.3) %>%
  pull(Protein) -> proteome_heps_30

d_binuc %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_binuc)) %>%
  filter(completeness > 0.5) %>%
  pull(Protein) -> proteome_heps_50

d_binuc %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_binuc)) %>%
  filter(completeness > 0.7) %>%
  pull(Protein) -> proteome_heps_70

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide

## QC of mono versus binucleated
# Number of proteins
d_binuc %>%
  group_by(cell_ID, Classification) %>%
  drop_na(int_core) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Classification, y = n, fill = Classification))+
  geom_boxplot() +
  scale_fill_manual(values = viridis(4)[2:3])+
  theme_classic() +
  labs(y = "Number of proteins") -> plot_binuc_nprot

# Intensity of core proteome
d_binuc %>%
  group_by(cell_ID, Classification) %>%
  drop_na(int_core) %>%
  filter(Protein %in% proteome_complete_binuc) %>%
  ggplot(aes(x = Classification, y = int_core))+
  geom_boxplot() +
  labs(title = "Intensities of core proteome")

d_binuc %>%
  group_by(cell_ID, Classification) %>%
  drop_na(int_core) %>%
  filter(Protein %in% proteome_complete_binuc) %>%
  group_by(Classification) %>%
  summarise(median = median(int_core))

# Intensities mono/binuc
d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  ggplot(aes(x = Classification, y = int_core)) +
  geom_boxplot() +
  labs(title = "Intensities of 30% complete proteome")

d_binuc %>%
  group_by(cell_ID, Classification) %>%
  drop_na(int_core) %>%
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Classification) %>%
  summarise(median = median(int_core))

## Fluorescence
img_fluovalues %>%
  left_join(meta_binuc) %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  filter(Channel == "Alexa488") %>%
  ggplot(aes(x = Classification, y = Area / 1000, fill = Classification)) +
  geom_boxplot()+
  scale_y_continuous(limits = c(500, 1000))+
  scale_fill_manual(values = viridis(4)[2:3])+
  theme_classic() -> plot_binuc_area

ggsave("../output/Figures/Binucleation_area.pdf")

img_fluovalues %>%
  left_join(meta_binuc) %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  filter(Channel != "Alexa488") %>%
  ggplot(aes(fill = Classification, y = Median, x = Channel)) +
  geom_boxplot()+
  scale_y_continuous(limits = c(0, 900))+
  scale_fill_manual(values = viridis(4)[2:3])+
  theme_classic() -> plot_binuc_fluo

## PCA
meta_distances_bins_4 %>%
  column_to_rownames("cell_ID") -> meta_pca

d_binuc %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  filter(cell_ID %in% rownames(meta_pca)) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("Protein") -> d_complete_binuc

# Plotting functions

p_heps <- PCAtools::pca(d_complete_binuc, metadata = meta_pca[colnames(d_complete_binuc),], removeVar = 0.1)

PCAtools::biplot(p_heps ,
                 colby = 'Classification',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = F,
                 legendPosition = 'right')+
  scale_color_manual(values = viridis(4)[2:3])  -> plot_pca12

ggsave(plot_pca12, file = "../output/Figures/PCA_binucleation_PC12.pdf")

PCAtools::eigencorplot(p_heps,
             metavars = c("Classification", "ratio"))

PCAtools::biplot(p_heps ,
                 x = 'PC4', y = 'PC5',
                 colby = 'Classification',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = T,
                 encircleFill = F,
                 showLoadings = T,
                 legendPosition = 'right')+
  scale_color_manual(values = viridis(4)[2:3]) -> plot_pca46

ggsave(plot_pca46, file = "../output/Figures/PCA_binucleation_PC46.pdf")

## Statistics with rstatix

# Complete proteome
d_binuc %>%
  filter(Protein %in% proteome_complete_binuc) %>%
  group_by(Protein) %>%
  group_by(Protein) %>%
  t_test(int_core ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_complete

ggplot(wilcox_complete, aes(x = estimate, y = -log10(p.adj)))+
  geom_point()+
  labs(title = "Binuc vs mono of core proteome")

# > 30% proteome
d_binuc %>%
  #filter(ratio < 0.5) %>%
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Protein) %>%
  filter(bio_ID == "m5C") %>%
  t_test(int_core ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_heps_30

ggplot(wilcox_heps_30, aes(x = estimate, y = -log10(p.adj)))+
  geom_point() +
  labs(title = "Binuc vs mono of 30% complete proteome, no space")

scaling_factor <- wilcox_heps_30 %>%
  summarise(scale = mean(estimate)) %>%
  pull(scale)

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  mutate(scaled_int = ifelse(Classification == "Binuc", int_core + scaling_factor, int_core)) %>%
  group_by(Protein) %>%
  t_test(scaled_int ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_heps_30_scaled

ggplot(wilcox_heps_30_scaled, aes(x = -estimate, y = -log10(p.adj)))+
  geom_point()+
  labs(title = "Binuc vs mono of 30% complete proteome, no space, scaled")

wilcox_heps_30_scaled %>%
  summarise(scale = mean(estimate)) %>%
  pull(scale)

WebGestaltR_gsea_binuc_30_rescaled <- WebGestaltR(enrichMethod="GSEA", organism="mmusculus",
                                                  enrichDatabase = c("geneontology_Biological_Process","geneontology_Cellular_Component","geneontology_Molecular_Function"),
                                                  interestGene = wilcox_heps_30_scaled %>%
                                                    dplyr::select(Protein, estimate),
                                                  interestGeneType ="uniprotswissprot")

# > 50% proteome close to PV

d_binuc %>%
  filter(ratio < 0.5) %>%
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Protein) %>%
  group_by(Protein) %>%
  t_test(int_core ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_heps_50_PV

ggplot(wilcox_heps_50_PV, aes(x = estimate, y = -log10(p.adj)))+
  geom_point()

scaling_factor <- wilcox_heps_50_PV %>%
  summarise(scale = mean(estimate)) %>%
  pull(scale)

d_binuc %>%
  filter(Protein %in% proteome_heps_50) %>%
  filter(ratio < 0.3) %>%
  mutate(scaled_int = ifelse(Classification == "Binuc", int_core + scaling_factor, int_core)) %>%
  group_by(Protein) %>%
  t_test(scaled_int ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_heps_50_PV_scaled

ggplot(wilcox_heps_50_PV_scaled, aes(x = -estimate, y = -log10(p.adj)))+
  geom_point()+
  theme_classic()

wilcox_heps_50_PV_scaled %>%
  summarise(scale = mean(estimate)) %>%
  pull(scale)

# > 50% proteome close to CV

d_binuc %>%
  filter(ratio > 0.7) %>%
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Protein) %>%
  group_by(Protein) %>%
  t_test(int_core ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_heps_30_CV

ggplot(wilcox_heps_30_CV, aes(x = estimate, y = -log10(p.adj)))+
  geom_point()

scaling_factor <-  wilcox_heps_30_CV %>%
  summarise(scale = mean(estimate)) %>%
  pull(scale)

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  filter(bio_ID == "m5C") %>%
  filter(ratio > 0.7) %>%
  mutate(scaled_int = ifelse(Classification == "Binuc", int_core + scaling_factor, int_core)) %>%
  group_by(Protein) %>%
  t_test(scaled_int ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_heps_30_CV_scaled

ggplot(wilcox_heps_30_CV_scaled, aes(x = -estimate, y = -log10(p.adj)))+
  geom_point()+
  theme_classic()

wilcox_heps_30_CV_scaled%>%
  summarise(scale = mean(estimate)) %>%
  pull(scale)

## Exclusive proteins

library(UpSetR)

d_binuc %>%
  distinct(Classification, Protein) %>%
  filter(Classification == "Mono") %>%
  pull(Protein) -> protein_mono

d_binuc %>%
  distinct(Classification, Protein) %>%
  filter(Classification == "Binuc") %>%
  pull(Protein) -> protein_binuc

listInput <- list(mono = protein_mono,
                  binuc = protein_binuc)

pdf("../output/Figures/Binucleation_Upset.pdf")
upset(fromList(listInput), order.by = "freq")
dev.off()

setdiff(protein_binuc,protein_mono) -> specific_binuc
setdiff(protein_mono,protein_binuc) -> specific_mono

meta_pg %>%
  mutate(mono_specific = Protein %in% specific_mono) %>%
  mutate(binuc_specific = Protein %in% specific_binuc) -> meta_pg_specific

## Plot protein intensities

d_binuc %>%
  filter(Protein %in% specific_mono) %>%
  group_by(Protein) %>%
  summarise(n = n()) %>%
  left_join(meta_pg) -> specific_mono_numbers

d_binuc %>%
  filter(Protein %in% specific_mono) %>%
  filter(Protein %in% (specific_mono_numbers %>% filter(n > 5) %>% pull(Protein))) %>%
  left_join(meta_pg) %>%
  ggplot(aes(x = Classification, y = int_core)) +
  geom_jitter(width = 0.1)+
  facet_wrap(.~Symbol)+
  theme_bw() -> plot_unique_mono

ggsave(plot_unique_mono, file = "../output/Figures/Binucleation_unique_mono.pdf")

d_binuc %>%
  filter(Protein %in% specific_binuc) %>%
  group_by(Protein) %>%
  summarise(n = n()) %>%
  left_join(meta_pg) -> specific_binuc_numbers

d_binuc %>%
  drop_na(Classification) %>%
  filter(Protein %in% specific_binuc) %>%
  filter(Protein %in% (specific_binuc_numbers %>% filter(n > 5) %>% pull(Protein))) %>%
  left_join(meta_pg) %>%
  ggplot(aes(x = Classification, y = int_core)) +
  geom_jitter(width = 0.1)+
  facet_wrap(.~Symbol)+
  theme_bw() -> plot_unique_binuc

ggsave(plot_unique_mono, file = "../output/Figures/Binucleation_unique_binuc.pdf")

## WebGestalt

WebGestaltR_ora_binuc <- WebGestaltR(enrichMethod="ORA", organism="mmusculus",
                                     enrichDatabase = c("geneontology_Biological_Process","geneontology_Cellular_Component","geneontology_Molecular_Function"),
                                     interestGene = specific_binuc,
                                     interestGeneType ="uniprotswissprot",
                                     referenceGene = meta_pg$Protein,
                                     referenceGeneType = "uniprotswissprot")

WebGestaltR_ora_mono <- WebGestaltR(enrichMethod="ORA", organism="mmusculus",
                                    enrichDatabase = c("geneontology_Biological_Process","geneontology_Cellular_Component","geneontology_Molecular_Function"),
                                    interestGene = specific_mono,
                                    interestGeneType ="uniprotswissprot",
                                    referenceGene = meta_pg$Protein,
                                    referenceGeneType = "uniprotswissprot")


# Pull cells that are positive for POI
d_binuc %>%
  left_join(meta_pg) %>%
  filter(Symbol %in% "Zc3h4") %>%
  pull(cell_ID) -> cells_cav1

# Perform statistics
meta_distances_bins_4 %>%
  mutate(COI = (cell_ID %in% cells_cav1)) %>%
  column_to_rownames("cell_ID") -> meta_tmp

d_binuc %>%
  #filter(Protein %in% proteome_heps_30) %>%
  filter(cell_ID %in% rownames(meta_tmp)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide_tmp

design <- model.matrix(~meta_tmp[colnames(d_wide_tmp),]$COI)

fit <- lmFit(d_wide_tmp, design)
fit <- eBayes(fit)
limma_cav1 <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  mutate(COI = (cell_ID %in% cells_cav1)) %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ COI,
              p.adjust.method = "fdr",
              detailed = T,
              ref.group = "FALSE") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox_cav1

d_binuc %>%
  filter(Protein == "P11930") %>%
  mutate(COI = (cell_ID %in% cells_cav1)) -> test
  ggplot(aes(x = COI, y = int_core))+
  geom_boxplot()
  

## Is it a Q value issue?
  
refquant <- read_tsv("Z:/Florian_Schober/scDVP/Reports/Liver-single-cell/DIANN/TO_CONSTANTIN__scDVP_full_report.tsv")

refquant %>%
  filter(Protein.Group %in% c("Q6ZPZ3","P16460")) %>%
  ggplot(aes(y = Channel.Q.Value, color = Protein.Group))+
  geom_density()


directlfq <- read_tsv("Z:/Florian_Schober/scDVP/Reports/Liver-single-cell/DIANN/TO_CONSTANTIN__scDVP_full_report.tsv.refquant_corrected.directLFQ.protein_intensities.tsv")

directlfq %>%
  gather(R.FileName, int, !protein) %>%
  mutate(int = ifelse(int == 0, NA, int)) %>%
  mutate(run_ID = paste(str_replace_all(str_replace_all(R.FileName, ".*scDVP_", ""), "_.*", ""),
                        str_replace_all(str_replace_all(R.FileName, ".*AID8_", ""), "_.*", ""),
                        sep = "_")) %>%
  mutate(label = str_replace(R.FileName, ".*_", "")) %>%
  mutate(label = str_replace_all(label, "\\(", ""), label = str_replace_all(label, "\\)", "")) %>%
  mutate(cell_ID = paste(run_ID, label, sep = "_")) %>%
  mutate(bio_ID = str_replace(cell_ID, "_.*", "")) %>%
  mutate(tube_ID = str_replace(run_ID, ".*_", "")) %>%
  dplyr::rename(Protein = protein) -> d_directlfq

d_directlfq %>%
  left_join(meta_distances_bins_4) %>%
  mutate(int_core = log2(as.numeric(as.character(int)))) -> d_binuc

d_binuc %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_binuc)) %>%
  filter(completeness > 0.3) %>%
  pull(Protein) -> proteome_heps_30

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  filter(cell_ID %in% SA_incl_binuc) %>%
  #filter(bio_ID == "m5C") %>%
  #filter(ratio > 0.7) %>%
  group_by(Protein) %>%
  t_test(int_core ~ Classification,
         ref.group = "Mono",
         detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_direclfq








img_fluovalues %>%
  mutate(Zc3h4_pos = cell_ID %in% cells_zc3h4) %>%
  filter(Channel == "CFP") %>%
  ggplot(aes(x = Zc3h4_pos, y = Area)) +
  geom_boxplot()

## Statistics
meta_distances_bins_4 %>%
  group_by(bin, Classification) %>%
  summarise(n = n())

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ Classification,
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox

d_binuc %>%
  
  filter(Protein %in% proteome_heps_30) %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ Classification,
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_bin1

## Limma statistics
design <- model.matrix(~meta_binuc[colnames(d_wide),]$Classification)

fit <- lmFit(d_wide, design)
fit <- eBayes(fit)
limma_binuc <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

ggplot(data = limma_binuc, aes(x = logFC, y = -log10(adj.P.Val)))+
  geom_point()

ggplot(data = limma_binuc, aes(y = logFC, x = 1))+
  geom_boxplot()

## PCA
d_binuc %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("Protein") -> d_complete_binuc

## Plotting functions

# - Loadings on PC1/PC2
p_heps <- PCAtools::pca(d_complete_binuc, metadata = meta_binuc[colnames(d_complete_binuc),], removeVar = 0.1)

PCAtools::biplot(p_heps ,
                 colby = 'Classification',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = F,
                 encircleFill = F,
                 showLoadings = T,
                 shape = 'Classification')







d_binuc %>%
  left_join(meta_pg) %>%
  filter(Symbol %in% (specific_mono_numbers %>% filter(n > 5) %>% pull(Symbol))) %>%
  group_by(cell_ID) %>%
  summarise(n = n()) %>%
  left_join(meta_distances_bins_4) -> specific_mono_cells

## Zc3h4
d_binuc %>%
  left_join(meta_pg) %>%
  filter(Symbol %in%  (specific_mono_numbers %>% filter(n > 5) %>% pull(Symbol))) %>%
  group_by(cell_ID) %>%
  left_join(meta_distances_bins_4) -> specific_mono_top3

# Pull cells that are positive for Zc3h4
d_binuc %>%
  left_join(meta_pg) %>%
  filter(Symbol %in% "Zc3h4") %>%
  pull(cell_ID) -> cells_zc3h4

# Perform statistics
meta_distances_bins_4 %>%
  #filter(bin %in% c(1,2,3)) %>%
  filter(Classification == "Binuc" | cell_ID %in% cells_zc3h4) %>%
  mutate(COI = (cell_ID %in% cells_zc3h4)) %>%
  column_to_rownames("cell_ID") -> meta_tmp

d_binuc %>%
  filter(Protein %in% proteome_heps_30) %>%
  filter(cell_ID %in% rownames(meta_tmp)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide_tmp

design <- model.matrix(~meta_tmp[colnames(d_wide_tmp),]$COI)

fit <- lmFit(d_wide_tmp, design)
fit <- eBayes(fit)
limma_tmp <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

# Bring in intensities



 # Wilcox
d_binuc %>%
  drop_na(int_core) %>%
  filter(Classification == "Mono") %>%
  mutate(COI = (cell_ID %in% cells_zc3h4)) %>%
  filter(Protein %in% proteome_complete_binuc) %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ COI,
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") ->  wilcox_pgrmc2

d_binuc %>%
  mutate(int_core = replace(int_core, is.na(int_core), 0)) %>%
  filter(Protein %in% specific_binuc) %>%
  left_join(meta_pg) %>%
  ggplot(aes(x = Classification, y = int_core)) +
  geom_jitter(width = 0.1)+
  facet_wrap(.~Symbol)

d_binuc %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc %>% rownames_to_column("cell_ID")) %>%
  filter(Protein %in% specific_binuc) %>%
  ggplot(aes(x = Classification, y = int_core, color = cell_ID)) +
  geom_jitter(width = 0.1)+
  facet_wrap(.~Protein)

meta_distances_bins %>%
  rownames_to_column("cell_ID") -> meta_distances

d_binuc %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc %>% rownames_to_column("cell_ID")) %>%
  filter(Protein %in% specific_binuc) %>%
  group_by(cell_ID) %>%
  summarise(n = n()) -> cell_hunting_binuc

d_binuc %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc %>% rownames_to_column("cell_ID")) %>%
  filter(Protein %in% specific_binuc) %>%
  filter(cell_ID == "m5C_26_target4") %>%
  left_join(meta_pg)

d_binuc %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc %>% rownames_to_column("cell_ID")) %>%
  filter(Protein %in% specific_mono) %>%
  group_by(cell_ID) %>%
  summarise(n = n()) -> cell_hunting_mono

## Spatial component
meta_binuc %>%
  rownames_to_column("cell_ID") %>%
  left_join(meta_distances) %>%
  mutate(bin = round((bin+1)/2, digits = 0)) %>%
  drop_na(bin) -> meta_binuc_distances

for(i in c(min(meta_binuc_distances$bin):max(meta_binuc_distances$bin))){
  meta_binuc_distances %>%
    filter(bin == i) %>%
    column_to_rownames("cell_ID") -> meta_tmp
  
  d_binuc %>%
    filter(Protein %in% proteome_heps_50) %>%
    filter(cell_ID %in% rownames(meta_tmp)) %>%
    dplyr::select(cell_ID, int_core, Protein) %>%
    spread(cell_ID, int_core) %>%
    arrange(Protein) %>%
    column_to_rownames("Protein") -> d_wide_tmp
  
  design <- model.matrix(~meta_tmp[colnames(d_wide_tmp),]$Classification)
  
  fit <- lmFit(d_wide_tmp, design)
  fit <- eBayes(fit)
  limma_tmp <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
    rownames_to_column("Protein") %>%
    left_join(meta_pg) %>%
    filter(adj.P.Val < 0.05) %>%
    pull(Protein)
  
  assign(paste("limma", i, sep = "_"), limma_tmp)
}

d_binuc %>%
  left_join(meta_binuc_distances) %>%
  filter(Protein == "O08807") %>%
  ggplot(aes(x = Classification, y = int_core, color = as.character(bin)))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), width = 0.1)+
  facet_wrap(.~bin)

meta_pg %>%
  filter(Protein %in% limma_2)


WebGestaltR_gsea_tmp <- WebGestaltR(enrichMethod="GSEA", organism="mmusculus",
                                    enrichDatabase = c("geneontology_Biological_Process","geneontology_Cellular_Component","geneontology_Molecular_Function"),
                                    interestGene = data.frame(limma_binuc %>% dplyr::select(ENSEMBL, logFC) %>% arrange(logFC)),
                                    interestGeneType ="ensembl_gene_id")


###########################
## Without normalization ##
###########################

binuc_prior_norm <- d_all %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  left_join(meta_binuc) %>%
  mutate(int = log2(int))

ggplot(binuc_prior_norm, aes(x = Classification, y = log2(int)))+
  geom_boxplot()

binuc_prior_norm %>%
  mutate(int = log2(int)) %>%
  group_by(Classification) %>%
  summarise(median = median(int, na.rm = T))

binuc_prior_norm %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  drop_na(int) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(rownames(meta_binuc))) %>%
  filter(completeness > 0.7) %>%
  pull(Protein) -> proteome_heps_50

binuc_prior_norm %>%
  filter(Protein %in% proteome_heps_50) %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  dplyr::select(cell_ID, int, Protein) %>%
  spread(cell_ID, int) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide

## Intensities mono/binuc
binuc_prior_norm %>%
  filter(Protein %in% proteome_heps_50) %>%
  dplyr::select(cell_ID, int, Protein) %>%
  left_join(meta_binuc %>% rownames_to_column("cell_ID")) %>%
  ggplot(aes(x = Classification, y = int)) +
  geom_boxplot()

## Limma statistics
design <- model.matrix(~meta_binuc[colnames(d_wide),]$Classification)

fit <- lmFit(d_wide, design)
fit <- eBayes(fit)
limma_binuc <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg)

ggplot(data = limma_binuc, aes(x = logFC, y = -log10(adj.P.Val)))+
  geom_point()

ggplot(data = limma_binuc, aes(y = logFC, x = 1))+
  geom_boxplot()


## Identify via number of missing values
d_binuc %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide


d_wide %>%
  rownames_to_column("Protein") %>%
  gather(cell_ID, value, !Protein) %>%
  mutate(value = ifelse(is.na(value), 0, 1)) %>%
  left_join(meta_binuc_distances) %>%
  group_by(Protein, Classification) %>%
  summarise(proportion = sum(value)/n()) %>%
  spread(Classification, proportion) %>%
  mutate(diff = `Binuc`- `Mono`) %>%
  left_join(meta_pg) -> test

## Rstatix
library(rstatix)

d_binuc %>%
  filter(Protein %in% proteome_heps_50) %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc_distances) %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ Classification,
              ref.group = "Mono",
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox

d_binuc %>%
  filter(Protein == "Q4LDG0") %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc_distances) %>%
  ggplot(aes(x = Classification, y = int_core))+
  geom_boxplot()

d_binuc %>%
  filter(Protein == "Q4LDG0") %>%
  filter(cell_ID %in% rownames(meta_binuc)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_binuc_distances) %>%
  ggplot(aes(x = as.factor(bin), y = int_core, color = Classification))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.2))

meta_binuc_distances %>%
  group_by(Classification, bin) %>%
  summarise(n = n())
  ggplot(aes(x = bin, y = ))
  
  d_binuc %>%
    filter(Protein %in% proteome_heps_50) %>%
    filter(cell_ID %in% rownames(meta_binuc)) %>%
    dplyr::select(cell_ID, int_core, Protein) %>%
    left_join(meta_binuc_distances) %>%
    group_by(Protein, bin) %>%
    wilcox_test(int_core ~ Classification,
                ref.group = "Mono",
                p.adjust.method = "fdr",
                detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") -> wilcox
  
## Try different number of classes

meta_distances_bins_4 %>%
  group_by(bin, Classification) %>%
  summarise(n = n())

binOI = 1
SA_incl_bin = meta_distances_bins_4 %>%
  filter(bin == binOI) %>%
  pull(cell_ID)

d_binuc %>%
  filter(bin == binOI) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/(length(SA_incl_bin))) %>%
  right_join(d_binuc) %>%
  filter(bin == binOI) %>%
  filter(completeness > 0.7) %>%
  group_by(Protein, bin) %>%
  wilcox_test(int_core ~ Classification,
              ref.group = "Mono",
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox_bin1

binOI = 2
SA_incl_bin = meta_distances_bins_4 %>%
  filter(bin == binOI) %>%
  pull(cell_ID)

d_binuc %>%
  filter(bin == binOI) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/(length(SA_incl_bin))) %>%
  right_join(d_binuc) %>%
  filter(bin == binOI) %>%
  filter(completeness > 0.7) %>%
  group_by(Protein, bin) %>%
  wilcox_test(int_core ~ Classification,
              ref.group = "Mono",
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox_bin2

binOI = 3
SA_incl_bin = meta_distances_bins_4 %>%
  filter(bin == binOI) %>%
  pull(cell_ID)

d_binuc %>%
  filter(bin == binOI) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/(length(SA_incl_bin))) %>%
  right_join(d_binuc) %>%
  filter(bin == binOI) %>%
  filter(completeness > 0.7) %>%
  group_by(Protein, bin) %>%
  t_test(int_core ~ Classification,
              ref.group = "Mono",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  left_join(meta_pg) -> wilcox_bin3

ggplot(data = wilcox_bin3, aes(x = estimate, y = -log10(p.adj)))+
  geom_point()

binOI = 4
SA_incl_bin = meta_distances_bins_4 %>%
  filter(bin == binOI) %>%
  pull(cell_ID)

d_binuc %>%
  filter(bin == binOI) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/(length(SA_incl_bin))) %>%
  right_join(d_binuc) %>%
  filter(bin == binOI) %>%
  filter(completeness > 0.7) %>%
  group_by(Protein, bin) %>%
  wilcox_test(int_core ~ Classification,
              p.adjust.method = "fdr",
              detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox_bin4


## Normalization with sctransform
library(sctransform)

d_wide[is.na(d_wide)] = 0

normalized_data <- sctransform::vst(as.matrix(d_wide))$y
as.data.frame(normalized_data) %>%
  rownames_to_column("Protein") %>%
  gather(cell_ID, int_core, !Protein) %>%
  left_join(meta_distances_bins_4) %>%
  filter(cell_ID %in% SA_incl_binuc) -> sc_long

sc_long %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ Classification,
              p.adjust.method = "fdr",
              detailed = T,
              ref.group = "Mono") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox_sc

ggplot(data = wilcox_sc, aes(x = estimate, y = -log10(p)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))

## Comparing to limma
d_wide %>%
  rownames_to_column("Protein") %>%
  gather(cell_ID, int_core, !Protein) %>%
  left_join(meta_distances_bins_4) %>%
  filter(cell_ID %in% SA_incl_binuc) -> sc_long

sc_long %>%
  group_by(Protein) %>%
  wilcox_test(int_core ~ Classification,
              p.adjust.method = "fdr",
              detailed = T,
              ref.group = "Mono") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") -> wilcox_sc

ggplot(data = wilcox_sc, aes(x = estimate, y = -log10(p)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))


