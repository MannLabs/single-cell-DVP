###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())


## Read relevant data
load("../output/variables/statTable_all_n.R")
load("../output/variables/SA_incl_all.R")
load("../output/Variables/img_fluovalues.R")

order_of_measurement <- data.frame(order_ID = as.character(c(3,2,1)), bio_ID = c("m3B", "m4A", "m5C"))

## Plotting functions
statTable_all_n %>%
  left_join(order_of_measurement) %>%
  mutate(included = cell_ID %in% SA_incl_all) %>%
  left_join(img_fluovalues) %>%
  mutate(Area = Area) %>%
  ggplot(aes(x = Area, y = n, color = order_ID, pch = included)) +
  geom_point() +
  facet_wrap(.~order_ID) +
  theme_bw()+
  geom_smooth(formula = y ~ log(x), se = F, method = "lm") +
  labs(x = "Microdissection area (kpx2)", y = "#Protein IDs") +
  scale_color_manual(values = viridis(5)[2:4]) +
  theme_classic() +
  scale_shape_manual(values = c(13, 19)) -> plot_area_proteinIDs

## Save figure to file
ggsave(plot_area_proteinIDs, file = "../output/Figures/Protein-IDs_VS_Area.pdf", width = 7, height = 5)
