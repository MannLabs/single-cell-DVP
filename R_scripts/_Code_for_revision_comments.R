###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/statTable_all_n.R")
load("../output/variables/statTable_all.R")
load("../output/variables/SA_incl_all.R")

## Read image to MS crossreference file
img_meta <- read_csv("../data/meta_img-proteome.csv") %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))

## Read image quantification file
img_fluovalues <- read_csv("../data/meta_img-quantification.csv") %>%
  dplyr::rename(bio_ID = Slide) %>%
  left_join(img_meta) %>%
  drop_na(cell_ID)

## Plotting functions
statTable_all_n %>%
  filter(cell_ID %in% SA_incl_all) %>%
  left_join(img_fluovalues) %>%
  mutate(Area = Area / 10000) %>%
  drop_na(Area) %>%
  distinct(cell_ID, .keep_all = T) -> area

classes = 4

data.frame(cell_ID = area$cell_ID, Area = area$Area) %>%
  mutate(bins = cut_interval(Area, n = classes))  -> bins

bins %>%
  distinct(bins) %>%
  arrange(bins) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(bins) -> area_bins
  
statTable_all_n %>%
  filter(cell_ID %in% SA_incl_all) %>%
  left_join(area_bins) %>%
  drop_na(bin) %>%
  ggplot(aes(x = bins, y = n, fill = bins)) +
  geom_boxplot() +
  labs(x = "Microdissection area (kpx2)", y = "#Protein IDs") +
  scale_fill_manual(values = viridis(6)[2:5]) +
  scale_y_continuous(limits = c(0,3000))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> plot_dissection_area_binned

ggsave(plot_dissection_area_binned, file = "../output/Figures/Review_Protein-IDs_VS_Area_binned.pdf", width = 7, height = 5)
