###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Load additional functions and environment
source("./Functions/normalize_core.R")
source("./Functions/read_prediction.R")

## Load relevant data
load("../Output/variables/meta_pg.R")
load("../Output/variables/d_classes.R")
load("../Output/variables/d_full_set.R")

## Read image meta data
img_meta <- read_csv("../data/meta_img-proteome.csv") %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))

img_fluovalues <- read_csv("../data/meta_img-quantification.csv") %>%
  dplyr::rename(bio_ID = Slide) %>%
  left_join(img_meta) %>%
  drop_na(cell_ID)

## Read common contaminats
contaminants <- read_tsv("../data/221109_common-contaminants_cRAP.tsv")

## Read and wrangle proteomics data of new mouse
d_m1A_priorNorm <- d_full_set %>%
  filter(bio_ID %in% c("m1A")) %>%
  filter(!Protein %in% contaminants$Entry) %>%
  left_join(img_meta) %>%
  mutate(heps = Index <= 302)

statTable_m1A <- d_m1A_priorNorm %>%
  left_join(meta_pg) %>%
  drop_na(int)

## Filtering criteria
statTable_m1A %>%
  filter(label != "Dimethyl-n-0") %>%
  group_by(run_ID, label, cell_ID, bio_ID) %>%
  summarise(n = n()) -> statTable_m1A_n

CIn_upper_all <- median(statTable_m1A_n$n) + 3* sd(statTable_m1A_n$n)
CIn_lower_all <- median(statTable_m1A_n$n) - 1.5* sd(statTable_m1A_n$n)

statTable_m1A_n %>%
  filter(n > CIn_lower_all, n < CIn_upper_all) %>%
  pull(cell_ID) -> SA_incl_m1A

## Normalization to core proteome
d_m1A_priorNorm %>% 
  filter(label != "Dimethyl-n-0") %>%
  filter(cell_ID %in% SA_incl_m1A) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_m1A

d_m1A_norm <- normalize_core(data = d_m1A_priorNorm %>% filter(cell_ID %in% SA_incl_m1A), proteome_core = proteome_complete_m1A, SA_incl = SA_incl_m1A) %>%
  left_join(meta_pg)

# Prediction data
prediction_m1A <- read_prediction("../data/shape_probability_m1A.csv", bio_ID = "m1A")

prediction_m1A %>%
  left_join(d_m1A_norm) %>%
  drop_na(R.FileName) %>%
  distinct(Protein, cell_ID, .keep_all = T) %>%
  drop_na(int_weighted, int_core)  -> prediction_m1A_empirical

cor(prediction_m1A_empirical$int_core, prediction_m1A_empirical$int_weighted)

ggplot(prediction_m1A_empirical, aes(x = log2(int_weighted), y = int_core))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  #scale_y_continuous(limits = c(11.1,14.8))+
  #scale_x_continuous(limits = c(11.1,14.8))+
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1, lty = "dotted") +
  scale_fill_viridis()+
  theme_classic()+
  labs(x = "log2(Intensity predicted)", y = "log2(Intensity measured)") -> plot_m1A_predaccuracy

ggsave(file = "../Output/Figures/m1A_predaccuracy.pdf", plot_m1A_predaccuracy, width = 5, height = 5)
