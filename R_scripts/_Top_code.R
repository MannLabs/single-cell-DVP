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
              "viridis",
              "ggrepel",
              "pheatmap")

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
refquant <- "../Data/89f5921594c92799946e67301cc9f12b_scDVP_refquant_corrected.formatted_for_iq.tsv.maxlfq_iq_protein_intensities.tsv"
save(refquant, file = "../output/Variables/refquant.R")

## -- Call files in order of appearance
source("./Data-wrangling.R")
source("./Histone-levels.R")
source("./PCA_Endothelial.R")
source("./PCA_Hepatocytes.R")
source("./PCA_vs_geometric-distance.R")
source("./Completeness_vs_Intensity.R")
source("./Rank_versus_Intensity.R")
source("./CVs.R")
source("./Labelling-efficiency.R")
source("./Protein-IDs_vs_Area.R")
source("./Protein-IDs_vs_Runs.R")
source("./Pseudo-FACS.R")
source("./Prediction_class_proteomes.R")
source("./Prediction_new_mouse.R")
source("./Heatmap_global_distances.R")
source("./Clustered_tables.R")
source("./Heatmap_markers.R")
source("./Limma_statistics.R")
source("./Spatial_expression_top10.R")
source("./Pathway_Urea_Peroxisome.R")
source("./Shapiro.R")
