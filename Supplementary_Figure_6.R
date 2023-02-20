##### Description ####################

# This script contains all analyses of Supplementary Figure 6:

# Supplementary Figure 6a: Abundance normalized to 1 across 9 bins in Halpern et al. 6 (left, marker expression guided bins), and this scDVP data (PC1-guided bins).
# Supplementary Figure 6b: Intensity correlation of all hits and zonation markers.
# Supplementary Figure 6c: Intensity correlation highlighting outliers. 
# Supplementary Figure 6d: Correlation coefficient for targets across all bins, with multiple testing adjusted p value
# Supplementary Figure 6e: Gne set enrichment analysis on Pearson correlation coefficient.


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
  
}


##### Set Working Directory ####################

# this function locates the current script in the folder structure. 
# input data and output files will be saved following this path

setwd(dirname(getActiveDocumentContext()$path))   

# save by:
# "output/",Technology ,"/",current_dataset,"/

## Variables 

current_dataset <- "Supplementary_Figure_6"




##### Preparation of Data ##################

#- around 50% of the expressed liver genes are non-randomly spatially zonated
#- expression of genes over more than four orders of magnitude, including genes with low expression

# DVP data were received from Florian Rosenberger (20220922), RNAseq were downloaded from https://www.nature.com/articles/nature21065#Sec13


# Clusters in the RNA dataset that should be merged
## here, no merging will be done (in contrast to the last approaches)
## The merginng is a mock for consistency in the script here:

RNA_merged_1 <- c("Layer 1")
RNA_merged_2 <- c("Layer 2")
RNA_merged_3 <- c("Layer 3")
RNA_merged_4 <- c("Layer 4")
RNA_merged_5 <- c("Layer 5")
RNA_merged_6 <- c("Layer 6")
RNA_merged_7 <- c("Layer 7")
RNA_merged_8 <- c("Layer 8")
RNA_merged_9 <- c("Layer 9")

# summarize these clusters in a Vector 
# with this, we can also add or remove a cluster in the point above
Merging_vector <- c("RNA_merged_1","RNA_merged_2", "RNA_merged_3", "RNA_merged_4", "RNA_merged_5","RNA_merged_6", "RNA_merged_7", "RNA_merged_8", "RNA_merged_9")




# Load data 

DVP <- data.table::fread(file = paste0("data/",current_dataset,"/221123_scDVP_Protoemics-to-RNAseq__PC1distance_cluster_9_new.csv"),
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


RNAseq <- data.table::fread(file = paste0("data/",current_dataset,"/Halpern_S3_Average_Expression.csv"),
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





# Reshape DVP dataframe 

DVP <- DVP %>%
  select(!c("Genename","ENSEMBL"))

DVP$cluster_anchorprotein <- gsub('1', 'cluster_anchor_1', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('2', 'cluster_anchor_2', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('3', 'cluster_anchor_3', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('4', 'cluster_anchor_4', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('5', 'cluster_anchor_5', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('6', 'cluster_anchor_6', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('7', 'cluster_anchor_7', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('8', 'cluster_anchor_8', DVP$cluster_anchorprotein)
DVP$cluster_anchorprotein <- gsub('9', 'cluster_anchor_9', DVP$cluster_anchorprotein)

# long to wide dataframe

DVP <- DVP %>% 
  spread(cluster_anchorprotein, int_mean, fill = NA, convert = FALSE)


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

#DVP_normalized <- 2^(DVP_normalized) 


#DVP_normalized$Summed.intensities <- rowSums(DVP)

DVP_normalized$Summed.intensities = apply(DVP_normalized, 1, sum)

DVP_normalized <- DVP_normalized/DVP_normalized[,"Summed.intensities"]

DVP_normalized <- DVP_normalized %>%
  select(-c(Summed.intensities))







# Reshape RNASeq dataframe 

# Exchange 0 with NaN

RNAseq[RNAseq == 0] <- NA  


# Filter for empty rows in dataframe

RNAseq <- RNAseq %>%
  column_to_rownames("Gene Symbol")

RNAseq$NA.count <- rowSums(is.na(RNAseq))

RNAseq <- RNAseq %>%
  filter(NA.count != 9)


# Merge Cluster Columns 
# NA values are ignored in the mean


# Dataframe for collecting results from loop
RNAseq_mean <-  RNAseq[,!names(RNAseq) %in% c(colnames(RNAseq))]


for (i in Merging_vector) {
  
  print(paste0("Starting with Cluster ", i))
  
  # Make mean of columns defined in loop
  intData <- RNAseq %>%
    select(all_of(get(i))) %>%
    mutate(Mean = rowMeans(., na.rm = T)) %>%
    select(Mean) 
  
  # Rename column and merge with target dataframe
  
  names(intData)[names(intData) == "Mean"] <- paste0(i)
  
  RNAseq_mean <- merge(RNAseq_mean, intData, by=0, all=TRUE) 
  
  RNAseq_mean <- RNAseq_mean %>%
    as.data.frame() %>%
    column_to_rownames("Row.names")
  
}


# reverse cluster order for RNAseq
# RNAseq_mean

RNAseq_mean <- RNAseq_mean %>%
  rename(
    RNA_merged_9 = RNA_merged_1,
    RNA_merged_8 = RNA_merged_2,
    RNA_merged_7 = RNA_merged_3,
    RNA_merged_6 = RNA_merged_4,
    RNA_merged_5 = RNA_merged_5,
    RNA_merged_4 = RNA_merged_6,
    RNA_merged_3 = RNA_merged_7,
    RNA_merged_2 = RNA_merged_8,
    RNA_merged_1 = RNA_merged_9) %>%
  select(order(colnames(.)))


# Normalize intensities
# Intensities are summed per transcript, then the ratio of each intensity is calculated

RNAseq_normalized <- RNAseq_mean

# rows with NA values are excluded
RNAseq_normalized$Summed.intensities <- rowSums(RNAseq_normalized)

RNAseq_normalized <- RNAseq_normalized/RNAseq_normalized[,"Summed.intensities"]

RNAseq_normalized <- RNAseq_normalized %>%
  select(-c(Summed.intensities)) %>%
  rownames_to_column("Symbol")


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


# fix RNA names with several genes

for (p in 1:nrow(RNAseq_normalized)) {
  
  #print(p)
  #break}
  
  Symbol.oi <- RNAseq_normalized[p, "Symbol"]
  
  if (grepl(";", Symbol.oi, fixed=TRUE)) {
    
    #print(paste0("Correction is done for ", Symbol.oi))
    
    symbol.vec <- str_split(Symbol.oi, ";") %>% unlist()
    
    # find overlap between row vector and DVP data
    
    overlap <- intersect(symbol.vec, DVP.vect)
    
    if (length(overlap) == 1) {
      print(paste0("Correction of gene names done for ", overlap))
      
      # replace gene name in dataframe
      RNAseq_normalized[p,"Symbol"] <- paste0(overlap)
    }
    
    if (length(overlap) > 1) {
      warning(paste0("1 Gene with several protein isoforms identified. It is kept in both datasets, but won't match between them."))
      break
    }
    
  }
} #end of for loop

RNAseq_normalized <- RNAseq_normalized %>%
  column_to_rownames("Symbol")




# Clear workspace
rm(list=setdiff(ls(), c("DVP_normalized", "RNAseq_normalized", "current_dataset")))







##### Supplementary Figure 6a ##################

# Define targets of interest

target.vec <- c("Glul", "Cyp2e1", "Ass1", "Asl", "Alb", "Cyp2f2")


# Plot data for RNAseq merged df 

# filter dataframe
RNAseq_normalized_filtered <- RNAseq_normalized %>%
  filter(rownames(.) %in% target.vec) %>%
  rownames_to_column("Symbol")

# wide to long

RNAseq_normalized_filtered  <- RNAseq_normalized_filtered %>%
  gather(Merged_cluster, Expression, -c("Symbol"))


# plot 

RNA.plt <- ggplot(RNAseq_normalized_filtered, aes(x=Merged_cluster, y=Expression, group=Symbol)) +
  geom_line(aes(color=Symbol),size = 1)+
  geom_point(aes(color=Symbol), size = 2.5) +
  xlab(paste0("Cluster")) + 
  ylab(paste0("Normalized Abundance")) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        panel.grid = element_line(colour = "transparent"),
        axis.text.x = element_text(angle=45,hjust=1)) +
  labs(title = paste0("scRNAseq"))


RNA.plt

ggsave(file = paste0("output/",current_dataset,"/RNA_Line_plot_marker_proteins_normalized_dataset.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/RNA_Line_plot_marker_proteins_normalized_dataset.pdf"),
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
        panel.grid = element_line(colour = "transparent"),
        axis.text.x = element_text(angle=45,hjust=1)) +
  labs(title = paste0("scDVP"))


DVP.plt

ggsave(file = paste0("output/",current_dataset,"/DVP_Line_plot_marker_proteins_normalized_dataset.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/DVP_Line_plot_marker_proteins_normalized_dataset.pdf"),
       width = 7,
       height = 6)


# Clear workspace

rm(list=setdiff(ls(), c("DVP_normalized", "RNAseq_normalized", "current_dataset")))






##### Supplementary Figure 6b and c ##################


# Unify ClusterNames between Datasets 

# Make long dataframes for each dataset

DVP.long <- DVP_normalized %>%
  rownames_to_column("Symbol") %>%
  gather(Cluster, Intensity_DVP, -c(Symbol)) 

RNA.long <- RNAseq_normalized %>%
  rownames_to_column("Symbol") %>%
  gather(Cluster, Intensity_RNA, -c(Symbol)) 

# Create common cluster names

datasets <- c("DVP.long", "RNA.long")

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
    
    # For Cluster 9
    if (grepl(9, q, fixed=TRUE)) {
      k.int[,"Cluster"] <- gsub(q, 'C9', k.int[,"Cluster"])
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
  
  
  if (k == "RNA.long") {
    RNA.common <- k.int
  }
  
  
}






# Generate Correlation Scatter Plots 

intData <- merge(RNA.common, DVP.common, by = 0) %>%
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
  geom_abline(slope= 1, intercept=0 ,colour= "grey", linetype = "dashed", size = 1)


correlation.plt

ggsave(file = paste0("output/",current_dataset,"/Correlation_Cluster_all_after_matching_correction.png"),
       width = 7,
       height = 6)

ggsave(file = paste0("output/",current_dataset,"/Correlation_Cluster_all_after_matching_correction.pdf"),
       width = 7,
       height = 6)



# same plot but with outliers highlighted

outliers.DVP <- intData %>%
  filter(Intensity_DVP >= 0.35)

outliers.RNA <- intData %>%
  filter(Intensity_RNA >= 0.4)

correlation.plt.highlight <- ggplot() + 
  geom_point(data = intData, aes(x=intData[,1], y=intData[,2], color = Cluster),  alpha = 0.2) +
  geom_point(data = outliers.DVP, x=outliers.DVP[,1], y=outliers.DVP[,2], color = "red", size = 2) +
  geom_point(data = outliers.RNA, x=outliers.RNA[,1], y=outliers.RNA[,2], color = "blue", size = 2) +
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
  
  geom_text_repel(data = outliers.RNA,
                  aes(x=outliers.RNA[,1], y=outliers.RNA[,2],
                      label = rownames(outliers.RNA)),
                  size = 1.5,
                  color = "black",
                  max.overlaps = Inf,
                  segment.size =0.5)

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
    
    res <- cor.test(intData2[,"Intensity_RNA"], intData2[,"Intensity_DVP"], 
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
       subtitle = paste0("RNAseq vs. DVP \nNumber of significant hits: ", length(significants$Gene)))


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



rm(list=setdiff(ls(), c("DVP_normalized", "RNAseq_normalized", "current_dataset", "uplot.df")))




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


