###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Label d0
d0_eff <- read_tsv("../data/Labelling_efficiency/d0/evidence.txt") %>%
  dplyr::select(Sequence, Modifications, Experiment, Intensity) %>%
  mutate(conc = str_replace_all(Experiment, "_.*", "")) %>%
  mutate(rep = str_replace_all(str_replace_all(Experiment, ".*DDA_", ""), "_.*", "")) %>%
  mutate(modified = grepl("Dimeth", Modifications))

d0_eff %>%
  group_by(Experiment, Sequence, modified) %>%
  summarise(sum_channel = sum(Intensity)) %>%
  ungroup() %>%
  right_join(d0_eff) %>%
  spread(modified, sum_channel) %>%
  mutate(`FALSE` = ifelse(is.na(`FALSE`), 0, `FALSE`)) %>%
  mutate(`TRUE` = ifelse(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(sum = `FALSE` + `TRUE`) %>%
  mutate(fraction = `TRUE` / sum) %>%
  group_by(rep, conc) %>%
  summarise(efficiency = mean(fraction)) %>%
  mutate(label = "d0")-> d0_eff_summary

## Label d4
d4_eff <- read_tsv("../data/Labelling_efficiency/d4/evidence.txt") %>%
  dplyr::select(Sequence, Modifications, Experiment, Intensity) %>%
  mutate(conc = str_replace_all(Experiment, "_.*", "")) %>%
  mutate(rep = str_replace_all(str_replace_all(Experiment, ".*DDA_", ""), "_.*", "")) %>%
  mutate(modified = grepl("Dimeth", Modifications))

d4_eff %>%
  group_by(Experiment, Sequence, modified) %>%
  summarise(sum_channel = sum(Intensity)) %>%
  ungroup() %>%
  right_join(d4_eff) %>%
  spread(modified, sum_channel) %>%
  mutate(`FALSE` = ifelse(is.na(`FALSE`), 0, `FALSE`)) %>%
  mutate(`TRUE` = ifelse(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(sum = `FALSE` + `TRUE`) %>%
  mutate(fraction = `TRUE` / sum) %>%
  group_by(rep, conc) %>%
  summarise(efficiency = mean(fraction)) %>%
  mutate(label = "d4")-> d4_eff_summary

## Label d8
d8_eff <- read_tsv("../data/Labelling_efficiency/d8/evidence.txt") %>%
  dplyr::rename(Experiment = `Raw file`) %>%
  dplyr::select(Sequence, Modifications, Experiment, Intensity) %>%
  mutate(conc = str_replace_all(Experiment, "_.*", "")) %>%
  mutate(rep = str_replace_all(str_replace_all(Experiment, ".*DDA_", ""), "_.*", "")) %>%
  mutate(modified = grepl("Dimeth", Modifications))

d8_eff %>%
  group_by(Experiment, Sequence, modified) %>%
  summarise(sum_channel = sum(Intensity)) %>%
  ungroup() %>%
  right_join(d8_eff) %>%
  spread(modified, sum_channel) %>%
  mutate(`FALSE` = ifelse(is.na(`FALSE`), 0, `FALSE`)) %>%
  mutate(`TRUE` = ifelse(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(sum = `FALSE` + `TRUE`) %>%
  mutate(fraction = `TRUE` / sum) %>%
  group_by(rep, conc) %>%
  summarise(efficiency = mean(fraction)) %>%
  mutate(label = "d8")-> d8_eff_summary

## Combine all three labels
rbind(d0_eff_summary, d4_eff_summary, d8_eff_summary) %>%
  group_by(label) %>%
  summarise(mean = mean(efficiency), sd = sd(efficiency, na.rm = T)) -> efficiency
  
  ggplot(efficiency, aes(x = label, y = mean, fill = label))+
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = mean + sd, ymin = mean)) +
  geom_text(aes(label = round(mean, 3)), hjust = 5, color = "white", size = 4, angle = 90) +
  theme_bw()+
  scale_fill_manual(values = rev(viridis(5)[4:2]))+
  scale_y_continuous(limits = c(0,1)) -> plot_labelling

ggsave(plot_labelling, file = "../output/Figures/Labelling_efficiency.pdf", width = 5, height = 5)
