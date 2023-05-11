###########################
#### scDVP Figure Code ####
###########################

#### -- Script caller -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## -- Package setup

packages <- c("rstudioapi",
              "tidyverse",
              "WebGestaltR",
              "limma",
              "org.Mm.eg.db",
              "org.Hs.eg.db",
              "viridis",
              "ggrepel",
              "pheatmap",
              "circlize",
              "ComplexHeatmap")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

## -- Set Working Directory
setwd(dirname(getActiveDocumentContext()$path))

# check if the folder exists
lapply(c("../output/Figures", "../output/Variables/", "../output/Tables"), function(paths){
  if (!dir.exists(paths)) {
    # create the folder if it doesn't exist
    dir.create(paths)
    cat("Folder created: ", paths, "\n")
  } else {
    cat("Folder already exists: ", paths, "\n")
  }
})

invisible(file.remove(list.files("../output/Figures", full.names = TRUE)))
invisible(file.remove(list.files("../output/Variables", full.names = TRUE)))
invisible(file.remove(list.files("../output/Tables", full.names = TRUE)))

## -- Set RefQuant data file
refquant <- "../data/protein/proteintable_scDVP.tsv"
save(refquant, file = "../output/Variables/refquant.R")

## -- Call files in order of appearance
source("./_data-wrangling.R")

    # Supplementary Figure 1 and 2
    source("./SuppFig1and2_Five_shapes.R")

# Figure 2
source("./Fig2_Histone-levels.R")

# Supplementary Figure 3
    source("./SuppFig3_Completeness_vs_Intensity.R")
    source("./SuppFig3_Rank_versus_Intensity.R")
    source("./SuppFig3_CVs.R")
    source("./SuppFig3_Labelling-efficiency.R")
    source("./SuppFig3_Protein-IDs_vs_Area.R")
    source("./SuppFig3_Protein-IDs_vs_Runs.R")

    # Supplementary Figure 4
    source("./SuppFig4_PCA_Endothelial.R")
    source("./SuppFig4_PCA_Hepatocytes.R")

# Figure 3
source("./Fig3_PCA_vs_geometric-distance.R")
source("./Fig3_Heatmap_global_distances.R")
source("./Fig3_Heatmap_markers.R")
source("./Fig3_Limma_statistics.R")
source("./Fig3_Spatial_expression_top10.R")
source("./Fig3_Pathway_Urea_Peroxisome.R")
source("./Fig3_GSEA.R")
source("./Fig3_Subcellular_localisation.R")

    #Supplemetary Figure 5
    source("./SuppFig5_PCA_reductive.R")
    
    # Supplementary Figure 6
    source("./SuppFig6_Shapiro.R")
    
    # Supplementary Figure 7
    source("./_clustered_tables.R")
    source("./SuppFig7_Comparison_to_RNAseq_9_Clusters.R")
    source("./SuppFig7_Comparison_to_FACS_8_Clusters.R")

# Figure 4
source("./Fig4_Pseudo-FACS.R")
source("./Fig4_Prediction_class_proteomes.R")
source("./Fig4_Prediction_new_mouse.R")
source("./Fig4_Prediction_m4A.R")

# Supplementary Figure 8
source("./SuppFig8_PCA_kmeans.R")







