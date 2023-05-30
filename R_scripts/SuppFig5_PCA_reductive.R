###########################
#### scDVP Figure Code ####
###########################

#### -- Supplementary Figure S5 -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/d.R")
load("../output/variables/meta_distances.R")

## Reductive PCA
d %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("Protein") -> d_complete_heps

meta_heps <- meta_distances %>%
  filter(cell_ID %in% colnames(d_complete_heps)) %>%
  arrange(ratio) %>%
  column_to_rownames("cell_ID")

## Calculate all divisors of length of meta_heps
num <- nrow(meta_heps)
divisors <- c()

# Find possible divisors
for (i in 1:num) {
  if (num %% i == 0) {
    divisors <- c(divisors, i)
  }
}

#divisors <- divisors[seq(1, length(divisors), 2)]

# Define a function that assigns repetitive elements
concatenate <- function(num, groups){
  y = rep(c(1:(num/groups)), each = groups)
  return(y)
}

## Add columns to metadata that allow concatenation
for(i in divisors){
  meta_heps %>%
    mutate(new_var = concatenate(num, i), !!paste0("1:", formatC(i, width = 3, flag = "0"), sep = "") := new_var) -> meta_heps
}

d %>%
  dplyr::select(Protein, int_core, cell_ID) %>%
  spread(cell_ID, int_core) %>%
  filter(complete.cases(.)) %>%
  gather(cell_ID, int, !Protein) %>%
  left_join(meta_heps %>% rownames_to_column("cell_ID")) %>%
  drop_na(ratio) %>%
  gather(concat, counter, grep("^1:", names(.), value = TRUE)) -> d_concat

d_concat %>%
  group_by(concat, counter, Protein) %>%
  summarise(int_concat = log2(median(2^int))) %>%
  mutate(sample = paste(concat, counter, sep = "_")) -> d_concat_summary

as.data.frame(d_concat_summary) %>%
  dplyr::select(sample, Protein, int_concat) %>%
  spread(sample, int_concat) %>%
  column_to_rownames("Protein") -> d_concat_summary_wide
  
d_concat_summary %>%
  distinct(sample) %>%
  column_to_rownames("sample") -> d_concat_meta

## Plotting functions
p_concat <- PCAtools::pca(d_concat_summary_wide[,rownames(d_concat_meta)], metadata = d_concat_meta, removeVar = 0.1)
# 
# PCAtools::biplot(p_concat,
#                  colby = 'concat',
#                  hline = 0, vline = 0,
#                  labSize = 3,
#                  lab = NA,
#                  encircle = F,
#                  encircleFill = F,
#                  showLoadings = F,
#                  legendPosition = 'right',
#                  alpha = 0.6)+
#   scale_color_manual(values = viridis(length(divisors))) +
#   theme_classic()-> plot_pca_loadings

## Plot drop in IQR depending on PC
p_concat[["rotated"]][,1:5] %>%
  rownames_to_column("sample") %>%
  left_join(d_concat_meta %>% rownames_to_column("sample")) %>%
  gather(component, value, grep("^PC", names(.), value = T)) %>%
  group_by(concat, component) %>%
  summarise(q1 = quantile(value, 0.25), q3 = quantile(value, 0.75), iqr = abs(q3 - q1)) %>%
  ggplot(aes(x = concat, y = iqr, color = component, group = component)) +
  geom_line()+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = viridis(5)) -> p_iqr

p_concat[["rotated"]][,1:3] %>%
  rownames_to_column("sample") %>%
  left_join(d_concat_meta %>% rownames_to_column("sample")) %>%
  gather(component, value, grep("^PC", names(.), value = T)) %>%
  ggplot(aes(x = concat, y = value, fill = component)) +
  geom_hline(yintercept = 0) +
  geom_boxplot()+
  scale_fill_manual(values = viridis(5)[2:4])+
  theme_classic()  -> p_box_pc

## Plot PCA depending on concatenation
p_concat[["rotated"]][,1:5] %>%
  rownames_to_column("sample") %>%
  left_join(d_concat_meta %>% rownames_to_column("sample")) %>%
  gather(component, value, grep("^PC", names(.), value = T)) %>%
  group_by(concat) %>%
  mutate(is_min = counter == 1, is_max = counter == max(counter)) %>%
  spread(component, value) -> p_reductive

# Calculate maximum distance lines
for(i in unique(p_reductive$concat)){
  
  p_concat[["rotated"]][,1:2] %>%
    rownames_to_column("sample") %>%
    left_join(d_concat_meta %>% rownames_to_column("sample")) %>%
    filter(concat == i) %>%
    dplyr::select(PC1, PC2) -> points
  
  dist_points <- as.matrix(dist(points))
  
  max_dist <- max(dist_points)
  max_pair <- which(dist_points == max_dist, arr.ind = TRUE)
  point1 <- points[max_pair[1], ]
  point2 <- points[max_pair[2], ]
  
  if(i == unique(p_reductive$concat)[1]){
    p_reductive_lines <- data.frame(rbind(point1, point2), i)
  } else {
    p_reductive_lines <- rbind(p_reductive_lines, data.frame(rbind(point1, point2), i))
  }
}

# Show in PCA
ggplot(p_reductive)+
  geom_hline(yintercept = 0, lty = "dotted")+
  geom_vline(xintercept = 0, lty = "dotted")+
  geom_point(aes(x = PC1, y = PC2, fill = concat, size = concat, alpha = concat), pch = 21, color = "black")+
  scale_fill_manual(values = viridis(length(divisors), option = "viridis")[length(divisors):1])+
  scale_alpha_manual(values = seq(0.3, 1, by = 0.5/length(divisors)))+
  theme_classic()+
  #geom_line(data = p_reductive_lines, aes(x = PC1, y = PC2, group = i), lty = "dashed")+
  #geom_point(data = p_reductive_lines, aes(x = PC1, y = PC2, fill = i), pch = 21, color = "black", size = 5, alpha = 1)+
  facet_wrap(.~concat) -> p_reductive

# Extract variance explained
for(i in unique(d_concat_summary$concat)){
  
  if(i == max(unique(d_concat_summary$concat))) break
  
  as.data.frame(d_concat_summary) %>%
    filter(concat == i) %>%
    dplyr::select(sample, Protein, int_concat) %>%
    spread(sample, int_concat) %>%
    column_to_rownames("Protein") -> d_concat_summary_wide
  
  d_concat_summary %>%
    filter(concat == i) %>%
    distinct(sample) %>%
    column_to_rownames("sample") -> d_concat_meta
  
  ## Plotting functions
  p_concat <- PCAtools::pca(d_concat_summary_wide[,rownames(d_concat_meta)], metadata = d_concat_meta, removeVar = 0.1)
  
  if(i == unique(d_concat_summary$concat)[1]){
    d_variance <- data.frame(variance = p_concat$variance, concat = i) %>%
      rownames_to_column("PC")
  } else{
    d_variance <- rbind(d_variance, data.frame(variance = p_concat$variance, concat = i) %>%
                          rownames_to_column("PC"))
  }
}

ggplot(data = d_variance %>% filter(PC %in% c("PC1", "PC2", "PC3", "PC4")), aes(x = concat, y = variance, group = PC, color = PC))+
  geom_line()+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = viridis(6)[2:5])+
  theme_classic() -> p_variance

ggsave(p_variance, file = "../output/Figures/Reductive_Variance.pdf", width = 5, height = 5)
ggsave(p_reductive, file = "../output/Figures/Reductive_PCA.pdf", width = 5, height = 5)
ggsave(p_iqr, file = "../output/Figures/Reductive_IQR.pdf", width = 5, height = 5)
ggsave(p_box_pc, file = "../output/Figures/Reductive_Boxplot.pdf", width = 7, height = 5)
