###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure XX -- ####
library(tidyverse)
library(ggrepel)

## Read relevant data
load("./output/statTable_all_n.R")
load("./output/SA_incl_all.R")

## Read image to MS crossreference file
img_meta <- read_csv("./metaImage/meta_img-proteome.csv") %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))

## Read image quantification file
img_fluovalues <- read_csv("./input/meta_image-quantification.csv") %>%
  dplyr::rename(bio_ID = Slide) %>%
  left_join(img_meta) %>%
  drop_na(cell_ID)

## Plotting functions
statTable_all_n %>%
  filter(cell_ID %in% SA_incl_all) %>%
  left_join(img_fluovalues) %>%
  mutate(Area = Area / 10000) %>%
  ggplot(aes(x = Area, y = n, color = bio_ID)) +
  geom_point() +
  facet_wrap(.~bio_ID) +
  theme_bw()+
  geom_smooth(formula = y ~ log(x), se = F, method = "lm") +
  labs(x = "Microdissection area (kpx2)", y = "#Protein IDs") +
  scale_color_manual(values = viridis(5)[2:4]) +
  theme_classic() -> plot_area_proteinIDs

## Save figure to file
ggsave(plot_area_proteinIDs, file = "./output/Figures/Protein-IDs_VS_Area.pdf", width = 7, height = 5)
