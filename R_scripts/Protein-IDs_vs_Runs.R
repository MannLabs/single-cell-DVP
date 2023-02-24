###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/statTable_all.R")
load("../output/variables/SA_incl_all.R")
load("../output/variables/CIn_lower_all.R")

## Plotting functions
median <- statTable_all %>%
  filter(label != "Dimethyl-n-0") %>%
  mutate(inclusion = cell_ID %in% SA_incl_all) %>%
  group_by(cell_ID, label, inclusion) %>%
  filter(inclusion == TRUE) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  summarise(median = median(n)) %>%
  pull(median)

statTable_all %>%
  mutate(inclusion = cell_ID %in% SA_incl_all) %>%
  group_by(cell_ID, label, inclusion, bio_ID, run_ID) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = run_ID, y = n, color = bio_ID, pch = inclusion))+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0,3000))+
  theme_bw()+
  scale_shape_manual(values = c(13, 19)) +
  geom_hline(yintercept = median, lty = "dashed", size = 1) +
  annotate("label", y = 3000, x = 30, label = paste("Median PG:", median)) +
  #geom_hline(yintercept = CIn_upper, lty = "dotted", size = 0.5) +
  geom_hline(yintercept = CIn_lower_all, lty = "dotted", size = 0.5) +
  labs(y = "# Protein IDs") +
  scale_color_manual(values = viridis(5)[2:4]) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  geom_hline(yintercept = 0) -> plot_NProt

ggsave(plot_NProt, file = "../output/Figures/Runs_vs_protein-IDs.pdf", width = 8, height = 5)
