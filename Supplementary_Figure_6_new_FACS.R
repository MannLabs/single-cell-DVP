##### Description ####################

# This script contains all analyses of Supplementary Figure 6:

# Integration of DVP and FACS sorted proteomics data

##### Prepare Workspace ####################

cat("\014")
gc()
rm(list=ls())
options(stringsAsFactors = FALSE)

##### Package Setup ####################

# set to "yes" if package setup is required
setup <- "no"

if (setup == "yes") {
  
  library(rstudioapi) 
  library(magrittr)
  library(tidyverse)
  library(dplyr)
  library(viridis)       
  library(ggpubr)
  
  library(formattable)
  library(htmltools)
  library(webshot)
  library(ggrepel)
  library(UpSetR)
  
  library(FactoMineR)
  library(factoextra)
  library(WebGestaltR)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  
  library(ggsci)
  library(gprofiler2)
  library(plotly)
  library(limma)
  library(statmod)
  
  library(svMisc)
  
  library(ggridges)
  
  library(readxl)
  
}


##### Set Working Directory ####################

# this function locates the current script in the folder structure. 
# input data and output files will be saved following this path

setwd(dirname(getActiveDocumentContext()$path))   

# save by:
# "output/",Technology ,"/",current_dataset,"/

## Variables 

current_dataset <- "Supplementary_Figure_6_new_FACS"



# Load data 


DVP  <- data.table::fread(file = paste0("output/tables/Proteome_to_FACS_9PCbins.tsv"),             
                                    sep = "\t",
                                    #nrows = 5000,
                                    header = TRUE,
                                    stringsAsFactors = FALSE,
                                    fill = TRUE,
                                    integer64 = "numeric",
                                    na.strings = "Filtered",
                                    data.table = FALSE,
                                    check.names = FALSE, # due to version issues
                                    verbose = FALSE,
                                    nThread = 10,
                                    showProgress = TRUE)

FACS <- data.table::fread(file = paste0("data/Supplementary_Figure_6/Moshe_S3_FACS_Proteomics.csv"),
                            sep = ",",
                            #nrows = 5000,
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            fill = TRUE,
                            integer64 = "numeric",
                            na.strings = "Filtered",
                            data.table = FALSE,
                            check.names = FALSE, # due to version issues
                            verbose = FALSE,
                            nThread = 10,
                            showProgress = TRUE)




##### 1. Reshape FACS dataframe ##################

# prepare FACS dataframe to look like RNAseq dataframe

FACS <- FACS %>%
 select(ends_with("LFQ")|c("Gene names", "Protein IDs"))

# remove log 2

all_samples <- as.character(colnames(FACS[, !colnames(FACS) %in% c("Gene names", "Protein IDs")] ))
FACS[,all_samples] <- 2^(FACS[,all_samples]) 

# make mean of replicates

mean.df <- data.frame(Columns = all_samples) %>%
  mutate(Cluster = str_split_fixed(.$Columns, "_", 3)[,2])

FACS.mean <- FACS %>%
  column_to_rownames("Protein IDs") %>%
  select(-c("Gene names")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Cluster = mean.df$Cluster[match(rownames(.), mean.df$Columns)]) %>%
  rownames_to_column("Sample IDs") %>%
  group_by(Cluster) %>%
  summarise(across(-`Sample IDs`, mean, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Cluster") %>%
  t() %>%
  as.data.frame()


# match back gene names

FACS.mean$`Gene names` <- FACS$`Gene names`[match(rownames(FACS.mean), FACS$`Protein IDs`)]

FACS.mean <- FACS.mean %>%
  rownames_to_column("Protein") 


colnames(FACS.mean)[colnames(FACS.mean) == 'variable'] <- 'bin'
colnames(FACS.mean)[colnames(FACS.mean) == 'value'] <- 'int'
colnames(FACS.mean)[colnames(FACS.mean) == 'Gene names'] <- 'Symbol'


FACS.mean$Symbol[FACS.mean$Symbol == ''] <- 'unknown_symbol'



##### 2a. Reshape and normalize DVP dataframe ##################

DVP <- DVP %>%
  select(!c("ENSEMBL"))

# remove logarithm

#DVP$int <- 10*(DVP$int)


# rename column

DVP$bin <- gsub('1', 'cluster_anchor_1', DVP$bin)
DVP$bin <- gsub('2', 'cluster_anchor_2', DVP$bin)
DVP$bin <- gsub('3', 'cluster_anchor_3', DVP$bin)
DVP$bin <- gsub('4', 'cluster_anchor_4', DVP$bin)
DVP$bin <- gsub('5', 'cluster_anchor_5', DVP$bin)
DVP$bin <- gsub('6', 'cluster_anchor_6', DVP$bin)
DVP$bin <- gsub('7', 'cluster_anchor_7', DVP$bin)
DVP$bin <- gsub('8', 'cluster_anchor_8', DVP$bin)
DVP$bin <- gsub('9', 'cluster_anchor_9', DVP$bin)

# long to wide dataframe

DVP <- DVP %>% 
  spread(bin, int, fill = NA, convert = FALSE)


# Check for row duplicates in gene names


DVP$Symbol[DVP$Symbol == 'NA'] <- 'unknown_symbol'

test_vector <- DVP$Symbol[duplicated(DVP$Symbol)] 

if (length(test_vector) > 0) {
  
  print("Adjustment to fix gene name duplicates is done.")
  
  
  for (z in test_vector) {
    
    # merge protein ID and gene name
    intData <- DVP %>%
      filter(Symbol == z) %>%
      mutate(Symbol = paste0(.$Symbol, "_",.$Protein ))
    
    # replace in dataframe
    
    DVP <- DVP %>%
      filter(Symbol != z)
    
    DVP <- rbind(DVP, intData)
    
    print(paste0("Gene symbol changed for ", z))
    
  }
  
}

DVP <- DVP %>%
  select(!c("Protein")) %>%
  column_to_rownames("Symbol")


# Normalize intensities
# Intensities are summed per protein, then the ratio of each intensity is calculated

DVP_normalized <- DVP

# reverse log2

DVP_normalized <- 2^(DVP_normalized) 


#DVP_normalized$Summed.intensities <- rowSums(DVP)

DVP_normalized$Summed.intensities = apply(DVP_normalized, 1, sum)

DVP_normalized <- DVP_normalized/DVP_normalized[,"Summed.intensities"]

DVP_normalized <- DVP_normalized %>%
  select(-c(Summed.intensities))




##### 2b. Reshape and normalize FACS dataframe ##################

# the DVP step from above is repeated

test_vector <- FACS.mean$Symbol[duplicated(FACS.mean$Symbol)] 


if (length(test_vector) > 0) {
  
  print("Adjustment to fix gene name duplicates is done.")
  
  
  for (z in test_vector) {
    
    # merge protein ID and gene name
    intData <- FACS.mean %>%
      filter(Symbol == z) %>%
      mutate(Symbol = paste0(.$Symbol, "_",.$Protein ))
    
    # replace in dataframe
    
    FACS.mean <- FACS.mean %>%
      filter(Symbol != z)
    
    FACS.mean <- rbind(FACS.mean, intData)
    
    print(paste0("Gene symbol changed for ", z))
    
  }
  
}

FACS.mean <- FACS.mean %>%
  select(!c("Protein")) %>%
  column_to_rownames("Symbol")


# Normalize intensities
# Intensities are summed per protein, then the ratio of each intensity is calculated

FACS_normalized <- FACS.mean

FACS_normalized$Summed.intensities = apply(FACS_normalized, 1, sum)

FACS_normalized <- FACS_normalized/FACS_normalized[,"Summed.intensities"]

FACS_normalized <- FACS_normalized %>%
  select(-c(Summed.intensities))



#### 3. Fix plurality of gene names ############################

# Fix plurality of gene names
# Some Genes have several names:
# Here we loop through each row and if there are several names, 
# we choose the one that is also in the DVP dataset

# generate vector with IDs that appear in DVP dataset
DVP.vect <- rownames(DVP_normalized) 


# remove unknown gene symbols in DVP data
DVP.vect <- DVP.vect[!DVP.vect %in% grep(paste0("unknown_symbol"), DVP.vect, value = T)]
# and remove these genes from DVP dataframe
DVP_normalized <- DVP_normalized %>%
  rownames_to_column("Symbol") %>%
  filter(., !grepl('unknown_symbol', Symbol)) %>%
  column_to_rownames("Symbol")

warning('Proteins with Symbol names "unknown_symbol" were removed.')




# generate vector with IDs that appear in FACS dataset
FACS.vect <- rownames(FACS_normalized) 

# remove unknown gene symbols in FACS data
FACS.vect <- FACS.vect[!FACS.vect %in% grep(paste0("unknown_symbol"), FACS.vect, value = T)]
# and remove these genes from FACS dataframe
FACS_normalized <- FACS_normalized %>%
  rownames_to_column("Symbol") %>%
  filter(., !grepl('unknown_symbol', Symbol)) 

warning('Proteins with Symbol names "unknown_symbol" were removed.')




# fix FACS gene names with several genes

for (p in 1:nrow(FACS_normalized)) {
  
  #print(p)
  #break}
  
  Symbol.oi <- FACS_normalized[p, "Symbol"]
  
  if (grepl(";", Symbol.oi, fixed=TRUE)) {
    
    #print(paste0("Correction is done for ", Symbol.oi))
    
    symbol.vec <- str_split(Symbol.oi, ";") %>% unlist()
    
    # find overlap between row vector and DVP data
    
    overlap <- intersect(symbol.vec, DVP.vect)
    
    if (length(overlap) == 1) {
      print(paste0("Correction of gene names done for ", overlap))
      
      # replace gene name in dataframe
      FACS_normalized[p,"Symbol"] <- paste0(overlap)
    }
    
    if (length(overlap) > 1) {
      warning(paste0("1 Gene with several protein isoforms identified. It is kept in both datasets, but won't match between them."))
      break
    }
    
  }
} #end of for loop


FACS_normalized <- FACS_normalized %>%
  column_to_rownames("Symbol")




# before integrating the data, we now exclude cluster 5 from the DVP data as this is the central proteomic sognal that was omitted in the FACS data 

DVP_normalized <- DVP_normalized %>%
  select(-c("cluster_anchor_5"))


# adjust column names

DVP.colnames <- data.frame(s.names = colnames(DVP_normalized)) %>%
  mutate(new_clusters = 1:nrow(.))

names(DVP_normalized) <- DVP.colnames$new_clusters[match(names(DVP_normalized), DVP.colnames$s.names)]

# Clear workspace
rm(list=setdiff(ls(), c("DVP_normalized", "FACS_normalized", "current_dataset")))




##### Supplementary Figure 6a ##################

# Define targets of interest

target.vec <- c("Glul", "Cyp2e1", "Ass1", "Asl", "Alb", "Cyp2f2")


# Plot data for RNAseq merged df 

# filter dataframe
FACS_normalized_filtered <- FACS_normalized %>%
  filter(rownames(.) %in% target.vec) %>%
  rownames_to_column("Symbol")

# wide to long

FACS_normalized_filtered  <- FACS_normalized_filtered %>%
  gather(Merged_cluster, Expression, -c("Symbol"))


# plot 

FACS.plt <- ggplot(FACS_normalized_filtered, aes(x=Merged_cluster, y=Expression, group=Symbol)) +
  geom_line(aes(color=Symbol),size = 1)+
  geom_point(aes(color=Symbol), size = 2.5) +
  xlab(paste0("Cluster")) + 
  ylab(paste0("Normalized Abundance")) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid = element_line(colour = "transparent")) +
  labs(title = paste0("FACS Proteome")) 


FACS.plt



ggsave(file = paste0("output/",current_dataset,"/FACS_Line_plot_marker_proteins_normalized_dataset.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/FACS_Line_plot_marker_proteins_normalized_dataset.pdf"),
       width = 7,
       height = 6)



#### Plot data for DVP ################

# filter dataframe
DVP_normalized_filtered <- DVP_normalized %>%
  filter(rownames(.) %in% target.vec) %>%
  rownames_to_column("Symbol")

# wide to long

DVP_normalized_filtered  <- DVP_normalized_filtered %>%
  gather(Merged_cluster, Expression, -c("Symbol"))



# plot 

DVP.plt <- ggplot(DVP_normalized_filtered, aes(x=Merged_cluster, y=Expression, group=Symbol)) +
  geom_line(aes(color=Symbol),size = 1)+
  geom_point(aes(color=Symbol), size = 2.5) +
  xlab(paste0("Cluster")) + 
  ylab(paste0("Normalized abundance")) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid = element_line(colour = "transparent")) +
  labs(title = paste0("scDVP"))


DVP.plt

ggsave(file = paste0("output/",current_dataset,"/DVP_Line_plot_marker_proteins_normalized_dataset.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/DVP_Line_plot_marker_proteins_normalized_dataset.pdf"),
       width = 7,
       height = 6)



# try combined plot

# for this, first match and combine dataframes

FACS_normalized_filtered$Technique <- paste0("FACS")
DVP_normalized_filtered$Technique <- paste0("DVP")

line.comb <- rbind(FACS_normalized_filtered, DVP_normalized_filtered)


library(plotly)


fig <- plot_ly(line.comb, x = ~Merged_cluster, y = ~Expression, z = ~Symbol, split = ~Technique, color = ~Symbol, type = 'scatter3d', mode = 'lines',
               line = list(width = 4, dash = c("solid", "dash")))

fig



FACS_normalized_filtered$Merged_cluster <- as.numeric(str_split_fixed(FACS_normalized_filtered$Merged_cluster, "_", 3)[,3])

DVP_normalized_filtered$Merged_cluster <- as.numeric(str_split_fixed(DVP_normalized_filtered$Merged_cluster, "_", 3)[,3])


scene = list(camera = list(eye = list(x = -1.25, y = -1.25, z = 1.25)))

fig <- plot_ly(type = 'scatter3d', mode = 'lines',
               line = list(width = 4)) %>%
  add_trace(data = DVP_normalized_filtered, x = ~Merged_cluster, y = ~Symbol, z = ~Expression, name = 'DVP', type = 'scatter3d', mode = 'lines', split = ~Symbol, color = ~Symbol,
            line = list(shape = 'linear', width= 4, dash = "dash")) %>%
  add_trace(data = FACS_normalized_filtered, x = ~Merged_cluster, y = ~Symbol, z = ~Expression, name = 'RNAseq', type = 'scatter3d', mode = 'lines', split = ~Symbol,color = ~Symbol,
            line = list(shape = 'linear', width= 4)) %>%
  layout(title = "DVP - RNAseq", scene = scene)


fig


# try paired plots with grid


comb.plt <- ggplot(line.comb, aes(x=Merged_cluster, y=Expression, group=Technique)) +
  geom_line(aes(color=Technique),size = 0.5)+
  geom_point(aes(color=Technique), size = 1.5) +
  xlab(paste0("Cluster")) + 
  ylab(paste0("Normalized abundance")) +
  scale_colour_viridis_d(begin = 0.8, end = 0) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid = element_line(colour = "transparent")) +
  facet_wrap(~Symbol, ncol = 3) +
  ggtitle(paste0("Marker panels of liver zonation - DVP vs. FACS Proteome"))


comb.plt

ggsave(file = paste0("output/",current_dataset,"/Line_plot_marker_proteins_normalized_combined.png"),
       width = 7,
       height = 5)

ggsave(file = paste0("output/",current_dataset,"/Line_plot_marker_proteins_normalized_combined.pdf"),
       width = 7,
       height = 5)


# Clear workspace

rm(list=setdiff(ls(), c("DVP_normalized", "FACS_normalized", "current_dataset")))



##### Supplementary Figure 6b and c ##################


# Unify ClusterNames between Datasets 

# Make long dataframes for each dataset

DVP.long <- DVP_normalized %>%
  rownames_to_column("Symbol") %>%
  gather(Cluster, Intensity_DVP, -c(Symbol)) 

FACS.long <- FACS_normalized %>%
  rownames_to_column("Symbol") %>%
  gather(Cluster, Intensity_FACS, -c(Symbol)) 

# Create common cluster names

datasets <- c("DVP.long", "FACS.long")

for (k in datasets) {
  
  k.int <- get(k)
  
  #print(k)}
  for (q in unique(get(k)[,"Cluster"])) {
    
    #print(q)}
    
    # For Cluster 1
    if (grepl(1, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C1', k.int[,"Cluster"])
    }
    
    # For Cluster 2
    if (grepl(2, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C2', k.int[,"Cluster"])
    }
    
    # For Cluster 3
    if (grepl(3, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C3', k.int[,"Cluster"])
    }
    
    # For Cluster 4
    if (grepl(4, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C4', k.int[,"Cluster"])
    }
    
    # For Cluster 5
    if (grepl(5, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C5', k.int[,"Cluster"])
    }
    
    # For Cluster 6
    if (grepl(6, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C6', k.int[,"Cluster"])
    }
    
    # For Cluster 7
    if (grepl(7, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C7', k.int[,"Cluster"])
    }
    
    # For Cluster 8
    if (grepl(8, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C8', k.int[,"Cluster"])
    }
    
    
  }
  
  
  # Generate unique identifier for each cluster
  
  k.int$ID <- paste0(k.int$Symbol, "_", k.int$Cluster)
  
  k.int <- k.int %>%
    column_to_rownames("ID") %>%
    select(-c("Symbol", "Cluster"))
  
  # rename dataframe
  
  if (k == "DVP.long") {
    DVP.common <- k.int
  }
  
  
  if (k == "FACS.long") {
    FACS.common <- k.int
  }
  
  
}



# Generate Correlation Scatter Plots 

intData <- merge(FACS.common, DVP.common, by = 0) %>%
  mutate(Symbol = paste0(str_split_fixed(.[,"Row.names"], "_", 2)[, 1]),
         Cluster = paste0(str_split_fixed(.[,"Row.names"], "_", 2)[, 2])) %>%
  column_to_rownames("Row.names")


# Generate extra dataframe with targets of interes

intData.targeted <- intData %>%
  filter(Symbol %in% c("Glul", "Ass1", "Oct", "Arg1", "Cps1",  "Asl", "Cyp2f2", "Pck1", "Cyp2e1")) 



# plot
correlation.plt <- ggplot() + 
  geom_point(data = intData, aes(x=intData[,1], y=intData[,2], color = Cluster),  alpha = 0.2) +
  geom_point(data = intData.targeted, x=intData.targeted[,1], y=intData.targeted[,2], color = "black", size = 2) +
  xlab(paste0(colnames(intData[1]))) + 
  ylab(paste0(colnames(intData[2]))) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid = element_line(colour = "transparent")) +
  geom_smooth(data = intData, aes(x=intData[,1], y=intData[,2]), method=lm, se=FALSE, color = "red", linetype = "dashed") + 
  geom_smooth(data = intData.targeted, aes(x=intData.targeted[,1], y=intData.targeted[,2]), method=lm, se=FALSE, color = "black", linetype = "dashed") + 
  labs(title = paste0("All clusters"),
       subtitle = paste0("Overall Pearson Correlation: ", signif(cor(intData[,1], intData[,2], use = "pairwise.complete.obs"),3)," (red)\n", "Marker Pearson Correlation: ", signif(cor(intData.targeted[,1], intData.targeted[,2], use = "pairwise.complete.obs"),3)," (black)")) +
  geom_text_repel(data = intData.targeted,
                  aes(x=intData.targeted[,1], y=intData.targeted[,2],
                      label = rownames(intData.targeted)),
                  size = 1.5,
                  color = "black")+
  annotate("text",x=0.5,y=0.25,label=(paste0("slope== ",coef(lm(intData[,2]~intData[,1]))[2])),parse=TRUE, color = "red")+
  annotate("text",x=0.5,y=0.4,label=(paste0("slope== ",coef(lm(intData.targeted[,2]~intData.targeted[,1]))[2])),parse=TRUE, color = "black") +
  annotate("text",x=0.5,y=0.5,label=(paste0("slope== 1")),parse=TRUE, color = "grey") +
  geom_abline(slope= 1, intercept=0 ,colour= "grey", linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(0,0.9)) +
  scale_y_continuous(limits = c(0,0.9)) 


correlation.plt

ggsave(file = paste0("output/",current_dataset,"/Correlation_Cluster_all_after_matching_correction.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/Correlation_Cluster_all_after_matching_correction.pdf"),
       width = 7,
       height = 6)



# same plot but with outliers highlighted

outliers.DVP <- intData %>%
  filter(Intensity_DVP >= 0.5)

outliers.FACS <- intData %>%
  filter(Intensity_FACS >= 0.4)

correlation.plt.highlight <- ggplot() + 
  geom_point(data = intData, aes(x=intData[,1], y=intData[,2], color = Cluster),  alpha = 0.2) +
  geom_point(data = outliers.DVP, x=outliers.DVP[,1], y=outliers.DVP[,2], color = "red", size = 2) +
  geom_point(data = outliers.FACS, x=outliers.FACS[,1], y=outliers.FACS[,2], color = "blue", size = 2) +
  xlab(paste0(colnames(intData[1]))) + 
  ylab(paste0(colnames(intData[2]))) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid = element_line(colour = "transparent")) +
  labs(title = paste0("Outliers by Intensity")) +
  geom_text_repel(data = outliers.DVP,
                  aes(x=outliers.DVP[,1], y=outliers.DVP[,2],
                      label = rownames(outliers.DVP)),
                  size = 1.5,
                  color = "black",
                  max.overlaps = Inf,
                  segment.size = 0.5)+
  
  geom_text_repel(data = outliers.FACS,
                  aes(x=outliers.FACS[,1], y=outliers.FACS[,2],
                      label = rownames(outliers.FACS)),
                  size = 1.5,
                  color = "black",
                  max.overlaps = Inf,
                  segment.size =0.5)+
  scale_x_continuous(limits = c(0,0.9)) +
  scale_y_continuous(limits = c(0,0.9)) 


correlation.plt.highlight


ggsave(file = paste0("output/",current_dataset,"/Correlation_Cluster_all_after_matching_correction_Outliers.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/Correlation_Cluster_all_after_matching_correction_Outliers.pdf"),
       width = 7,
       height = 6)




##### Supplementary Figure 6d ##################

# Pearson Correlation Testing 

# row wise test

# list of each row entry

entry_list <- unique(intData$Symbol)

uplot.df <- data.frame(p.val = NA, corr.coef = NA, Symbol = NA)

# Define runnung number for rows
i <- 1

for (k in entry_list) {
  
  #print(k)
  #break}
  intData2 <- intData %>%
    filter(Symbol == k)
  
  try({
    
    res <- cor.test(intData2[,"Intensity_FACS"], intData2[,"Intensity_DVP"], 
                    method = "pearson")
    
    # Extract the p.value
    p.val <- res$p.value
    
    # Extract the correlation coefficient
    corr.coef <- res$estimate[1]
    
    
    uplot.df[i,"p.val"] <- paste0(p.val)
    uplot.df[i,"corr.coef"] <- paste0(corr.coef)
    uplot.df[i,"Symbol"] <- paste0(k)
    
    i <- i + 1
    
    print(paste0("Done ", k))
    
  })
  
}

# multiple sample testing for p value

uplot.df$p.adjust <- p.adjust(uplot.df$p.val, method = "BH")


significants <- uplot.df %>%
  filter(p.adjust <= 0.05) %>%
  rownames_to_column("Gene")


marker.df <- uplot.df %>%
  filter(Symbol %in% c("Glul", "Ass1", "Oct", "Arg1", "Cps1",  "Asl", "Cyp2f2", "Pck1", "Cyp2e1"))


# define top positively correlating markers

posCor.df <- significants %>%
  arrange(., desc(as.numeric(corr.coef))) %>%
  mutate(rank = 1:nrow(.)) %>%
  filter(rank <= 10)


# define top negatively correlating markers

negCor.df <- significants %>%
  arrange(., desc(as.numeric(corr.coef))) %>%
  mutate(rank = 1:nrow(.)) %>%
  filter(corr.coef < 0 ) %>%
  filter(rank > (length(.$rank)-10))

# plot 

uplot <- ggplot() +
  geom_point(data=uplot.df, aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val))), size = 2, alpha = 0.2)+
  geom_point(data=significants, aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val))), size = 2, color = "#990000", alpha = 0.2) +
  geom_point(data=marker.df, aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val))), size = 2, color = "orange", alpha = 1) +
  geom_point(data=posCor.df, aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val))), size = 2, color = "#990000", alpha = 0.2) +
  geom_point(data=negCor.df, aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val))), size = 2, color = "#990000", alpha = 0.2) +
  theme_bw(base_size = 12) + 
  xlab (paste("Correlation Coefficient", sep = "")) +
  ylab ("-log10 (p.value)") +
  geom_text_repel(data = marker.df,
                  aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val)), label = Symbol),
                  size = 4,
                  color = "orange") +
  geom_text_repel(data = posCor.df,
                  aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val)), label = Symbol),
                  size = 4,
                  color = "#990000",
                  hjust=0,
                  segment.size=1,
                  nudge_x=-0.5, 
                  direction="y",
                  segment.alpha = 0.5,
                  segment.size = 1) + 
  geom_text_repel(data = negCor.df,
                  aes(x=as.numeric(corr.coef), y=-log10(as.numeric(p.val)), label = Symbol),
                  size = 4,
                  color = "#990000",
                  hjust=0,
                  segment.size=1,
                  nudge_x=0.5, 
                  direction="y",
                  segment.alpha = 0.5,
                  segment.size = 1)+
  labs(title = paste0("Pearson Correlation Testing"),
       subtitle = paste0("FACS vs. DVP \nNumber of significant hits: ", length(significants$Gene)))


uplot

ggsave(file = paste0("output/",current_dataset,"/Pearson_Correlation_Testing_Uplot.png"),
       width = 6,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/Pearson_Correlation_Testing_Uplot.pdf"),
       width = 6,
       height = 6)

# Export as Excel list

uplot.df <- uplot.df %>%
  mutate(corr.coef = as.numeric(corr.coef)) %>%
  arrange(., desc(corr.coef))

write.csv(uplot.df, file = paste0("output/",current_dataset,"/UPlot_Enrichment_", current_dataset, ".csv"))



rm(list=setdiff(ls(), c("DVP_normalized", "FACS_normalized", "current_dataset", "uplot.df")))




##### Supplementary Figure 6e ##################

# GSEA Term Enrichment of proteins ranked by Pearson correlation 

# Prepare Enrichment lists
enrichment.vector <- uplot.df %>%
  select(c("Symbol","corr.coef")) %>%
  as.data.frame(.)%>%
  arrange(., desc(as.numeric(corr.coef)))

# make all genes capital letters

enrichment.vector$Symbol = toupper(enrichment.vector$Symbol)

# enrichment with two lists

try({
  
  suppressWarnings(
    WebGestaltR::WebGestaltR(
      enrichMethod = "GSEA",
      organism = "hsapiens",
      enrichDatabase = c("pathway_KEGG", "pathway_Reactome", "geneontology_Biological_Process_noRedundant"), 
      interestGene = enrichment.vector,
      interestGeneType = "genesymbol",
      fdrThr = 0.05,
      isOutput = FALSE) %>%
      dplyr::rename() -> webgestalt_out
  )
  
  
  webgestalt_out <- webgestalt_out %>%
    select(c("description", "enrichmentScore", "normalizedEnrichmentScore", "pValue", "FDR", "size", "userId", "database"))
  
  write.csv(webgestalt_out, file = paste0("output/",current_dataset,"/GSEA_Enrichment_", current_dataset, ".csv"))
  
  
})




# Term of interest 

int.term <- "Glutathione conjugation"

gp_mod <- webgestalt_out

int.vector <- unlist(str_split(c(gp_mod[which(gp_mod$description == int.term),][,"userId"]), ";"))

gp_mod.filtered <- DVP_normalized %>%
  filter(toupper(rownames(.)) %in% int.vector)


# HeatMap

# Heatmap color

col_hm = colorRamp2(c(min(gp_mod.filtered, na.rm = T), (max(gp_mod.filtered, na.rm = T)/2), max(gp_mod.filtered, na.rm = T)), c("white","#feb24c", "#e31a1c"))

hm <- Heatmap(as.matrix(gp_mod.filtered),
              name = "Normalized Abundance",
              cluster_rows = T,
              cluster_columns = F,
              #right_annotation = col_annot,
              #bottom_annotation = col_annot,
              rect_gp = gpar(col = "grey87", lwd = 1),
              show_column_names = T,
              col = col_hm, 
              column_title = paste(int.term))


# png
png(file = paste0("output/",current_dataset,"/HeatMap_Positive_Correlation_", current_dataset,"_", int.term, ".png"), width = 7, height = 5, units = "in", res= 400)


draw(hm,
     heatmap_legend_side = "right")
dev.off()


pdf(file = paste0("output/",current_dataset,"/HeatMap_Positive_Correlation_", current_dataset,"_", int.term, ".pdf"), width = 7, height = 5)


draw(hm,
     heatmap_legend_side = "right")
dev.off()


