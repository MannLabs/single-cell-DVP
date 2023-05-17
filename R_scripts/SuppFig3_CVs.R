###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Additional function
source("./Functions/cv.R")

## Read relevant data
load("../output/variables/d.R")
load("../output/Variables/img_fluovalues.R")
load("../output/Variables/meta_distances.R")
load("../output/Variables/SA_incl_all.R")

## Define number of classes
classes = 4

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> meta_distances_bins

meta_distances_bins %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bins) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> meta_distances_bins

## Group by size and spatial bin
classes = 3

data.frame(cell_ID = img_fluovalues$cell_ID, Area = img_fluovalues$Area) %>%
  filter(cell_ID %in% SA_incl_all) %>%
  mutate(bins = cut_interval(Area, n = classes))  -> bins

bins %>%
  distinct(bins) %>%
  arrange(bins) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(bins) -> area_bins

d %>%
  filter(cell_ID %in% SA_incl_all) %>%
  filter(heps == TRUE) %>%
  left_join(meta_distances_bins %>% rownames_to_column("cell_ID")) %>%
  mutate(int = 2^int_core) %>%
  drop_na(bin) %>%
  dplyr::rename(spatial_bin = bin) %>%
  left_join(area_bins) %>%
  group_by(Protein, bin, spatial_bin) %>%
  mutate(cv = cv(int)) %>%
  mutate(spatial_bin = abs(spatial_bin - classes - 1)) %>%
  ggplot(aes(x = as.factor(spatial_bin), y = cv, fill = as.factor(bin)))+
  geom_boxplot() +
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120, 20))+
  theme_classic() +
  scale_fill_manual(values = viridis(5)[2:4]) -> plot_cvs_bins

d %>%
  filter(cell_ID %in% SA_incl_all) %>%
  filter(heps == TRUE) %>%
  left_join(meta_distances_bins %>% rownames_to_column("cell_ID")) %>%
  mutate(int = 2^int_core) %>%
  drop_na(bin) %>%
  dplyr::rename(spatial_bin = bin) %>%
  left_join(area_bins) %>%
  dplyr::rename(area_bin = bin) %>%
  mutate(spatial_bin = abs(spatial_bin - 9)) %>%
  group_by(area_bin, spatial_bin) %>%
  distinct(cell_ID) %>%
  summarise(n = n()) -> summary_stats

## Save figure to file
ggsave(plot_cvs_bins, file = "../output/Figures/CVs_bins.pdf", width = 6, height = 5)
  
