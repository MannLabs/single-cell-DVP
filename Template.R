##### Description ####################

# This script contains all analyses of (Supplementary) Figure xxxx:

# Figure xa: Short legend description

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
  
}


##### Set Working Directory ####################

# this function locates the current script in the folder structure. 
# input data and output files will be saved following this path

setwd(dirname(getActiveDocumentContext()$path))   

# save by:
# "output/",Technology ,"/",current_dataset,"/

