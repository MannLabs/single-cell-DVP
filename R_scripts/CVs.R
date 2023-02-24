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
load("../output/variables/d_all.R")
load("../output/variables/d_all_norm.R")

## Plotting functions
d_all %>%
  mutate(norm = FALSE) %>%
  dplyr::rename(int_core = int) %>%
  full_join(d_all_norm %>% mutate(norm = TRUE) %>% mutate(int_core = 2^int_core)) %>%
  group_by(Protein, norm) %>%
  mutate(cv = cv(int_core)) %>%
  ggplot(aes(x = norm, y = cv, fill = norm))+
  geom_boxplot() +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10))+
  theme_classic()+
  scale_fill_manual(values = viridis(4)[2:3]) -> plot_cvs

## Save figure to file
ggsave("../output/Figures/CVs.pdf", width = 4, height = 5)