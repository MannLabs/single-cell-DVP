###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 3A -- ####

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

divisors <- divisors[seq(1, length(divisors), 6)]

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

PCAtools::biplot(p_concat,
                 colby = 'concat',
                 hline = 0, vline = 0,
                 labSize = 3,
                 lab = NA,
                 encircle = T,
                 encircleFill = F,
                 showLoadings = F,
                 legendPosition = 'right')+
  scale_color_manual(values = viridis(length(divisors))) +
  theme_classic()#-> plot_pca_loadings
# 
# p_concat_values <- p_concat[["rotated"]][,1:2]
# 
# # Compute pairwise distances
# d <- as.matrix(dist(p_concat_values))
# 
# # Find the pair of points with maximum distance
# max_dist <- max(d)
# max_pair <- which(d == max_dist, arr.ind = TRUE)
# point1 <- points[max_pair[1], ]
# point2 <- points[max_pair[2], ]
# 
# 
# 
# 
# 
# 
# 
# ggsave(plot_pca_loadings, file = "../Output/Figures/PCA_Loadings.pdf", width = 7, height = 6)
# 
# ### Combine every nth cell
# 
# 
# 
# 
# mtcars_mean <- mtcars %>%
#   summarise(mpg_mean = mean(c(lag(mpg), mpg), na.rm = TRUE))
# 
# # View the result
# mtcars_mean
# 
