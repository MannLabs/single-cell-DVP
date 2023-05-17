Florian Rosenberger et al., 2023 - in revision

# Spatial single-cell mass spectrometry defines zonation of the hepatocyte proteome

## Abstract

Single-cell proteomics by mass spectrometry (MS) is emerging as a powerful and unbiased method for the characterization of biological heterogeneity. So far, it has been limited to cultured cells, whereas an expansion of the method to complex tissues would greatly enhance biological insights. Here we describe single-cell Deep Visual Proteomics (scDVP), a technology that integrates high-content imaging, laser microdissection and multiplexed MS. scDVP resolves the context-dependent, spatial proteome of murine hepatocytes at a current depth of 1,700 proteins from a slice of a cell. Half of the proteome was differentially regulated in a spatial manner, with protein levels changing dramatically in proximity to the central vein. We applied machine learning to proteome classes and images, which subsequently inferred the spatial proteome from imaging data alone. scDVP is applicable to healthy and diseased tissues and complements other spatial proteomics or spatial omics technologies.

## Table of contents

1. [Data repository](#Data-repository)
2. [Results](#Results)
3. [R Scripts](#R-Scripts)
4. [GitHub Notes](#GitHub-Notes)  

## Data repository

Processed mass spectrometry raw data and other [input](/input/) files have been saved in the following folders:

- [Data and metadata](/data/)

## Results

Figures and result dataframes are saved in the [output](/output/) folder. 

- [Results - Figures](/output/Figures/)
- [Results - Figures](/output/Tables/)

## R scripts

Run "_Top_code.R" to execute the entire R code.

### Data wrangling scripts
- [Header file](R_scripts/_Top_code.R)
- [Filtering criteria](R_scripts/Data-wrangling.R)
- [Prediction proteome calculation](R_scripts/Prediction_class_proteomes.R)
- [Clustered tables](R_scripts/Data-wrangling.R)

### Scripts for main figures
#### Figure 2, Depth of single shape proteomes and estimation of nuclear compartment
- [Figure 2a, Abundance range and transcription factors](R_scripts/Fig2_Rank_versus_Intensity.R)
- [Figure 2b, Intensity of top four histones](R_scripts/Fig2_Histone-levels.R)

#### Figure 3, Single shape proteomes are accurate descriptors of zonated hepatocytes
- [Figure 3a, PCA with distance metric overlay](R_scripts/Fig3_PCA_vs_geometric-distance.R)
- [Figure 3b, Distance metric versus PC1](R_scripts/Fig3_PCA_vs_geometric-distance.R)
- [Figure 3c, Global heatmap of 70% complete proteome](R_scripts/Fig3_Heatmap_global_distances.R)
- [Figure 3d, Expression by location](R_scripts/Fig3_Heatmap_markers.R)
- [Figure 3e, Relative expression of top-10 differential proteins](R_scripts/Fig3_Spatial_expression_top10.R)
- [Figure 3f, Spatial Gene Set Enrichment Analysis](R_scripts/Fig3_GSEA.R)
- [Figure 3g, Subcellular compartments](R_scripts/Fig3_GSEA.R)
- [Figure 3h, OXPHOS and mitochondrial fatty acid metabolism](R_scripts/Fig3_OXPHOS.R)
- [Figure 3i/j, Urea cycle and peroxisomal fatty acid metabolism](R_scripts/Fig3_Pathway_Urea_Peroxisome.R)

#### Figure 4, Combining imaging and proteome data for a machine-learned model
- [Figure 4a, Pseudo-FACS plot](R_scripts/Fig4_Pseudo-FACS.R)
- [Figure 4b, Proteome bin versus staining intensity](R_scripts/Pseudo-FACS.R)
- [Figure 4e, Proteome prediction of an unseen section](R_scripts/Prediction_new_mouse.R)
- [Figure 4f, Proteome prediction on m4A](R_scripts/Fig4_Prediction_m4A.R)

### Scripts for supplementary figures
#### Supplementary Figure S1 and S2: Five-shape proteomes
- [Figure S1 and S2](R_scripts/SuppFig1and2_Five_shapes.R)
- 
#### Supplementary Figure S3: Performance overview of single shape proteomes
- [Figure S3a, Labelling efficiency](R_scripts/SuppFig3_Labelling-efficiency.R)
- [Figure S3b, Area of hepatocytes]()
- [Figure S3c, Protein IDs per biological replicate](R_scripts/SuppFig3_Protein-IDs_vs_Runs.R)
- [Figure S3d, Microdissected area versus protein IDs](R_scripts/SuppFig3_Protein-IDs_vs_Area.R)
- [Figure S3e, Data completeness versus intensity](R_scripts/SuppFig3_Completeness_vs_Intensity.R)
- [Figure S3f, Normalization versus CVs](R_scripts/SuppFig3_CVs.R)

#### Supplementary Figure S4: Dimensionality reduction of single shape data
- [Figure S4a, PCA with Ass1 expression levels](R_scripts/SuppFig4_PCA_Hepatocytes.R)
- [Figure S4b, PCA with Cyp2e1 expression levels](R_scripts/SuppFig4_PCA_Hepatocytes.R)
- [Figure S4c, Distance metric versus PC12](R_scripts/SuppFig4_PCA_vs_geometric-distance.R)
- [Figure S4d, PCA with loadings](R_scripts/SuppFig4_PCA_Hepatocytes.R)
- [Figure S4e, PCA including endothelial cells](R_scripts/SuppFig4_PCA_Endothelial.R)
- 
#### Supplementary Figure S5: Concatentation of shapes
- [Figure S4](R_scripts/SuppFig5_PCA_reductive.R)

#### Supplementary Figure S6: Normality checks of single shape data
- [Figure S5a, Volcano plot](R_scripts/SuppFig6_Volcano_plot.R)
- [Figure S5b, Shapiro-Wilk-Test](R_scripts/Shapiro.R)
- [Figure S5d, Relative expression of bottom-10 differential proteins](R_scripts/Shapiro.R)


#### Supplementary Figure S7: Additional data on subcellular localisation
- [Figure S7, Linear models](R_scripts/Fig3_Subcellular_localisation.R)

#### Supplementart Figure S8
- [Figure S7a-c, Comparison to RNAseq](R_scripts/SuppFig7_Comparison_to_RNAseq_9_Clusters.R)
- [Figure S7d-g, Comparison to FACS-based Proteomics](R_scripts/SuppFig7_Comparison_to_FACS_8_Clusters.R)

#### Supplementart Figure S9
- [Figure S7d-g, Comparison to FACS-based Proteomics](SuppFig9_PCA_kmeans.R


## GitHub Notes

**Clone Repository**

Navigate to your local GitHub folder and enter:

`>git clone https://github.com/MannLabs/single-cell-DVP.git`

**Transferring local copy to GitHUb**

Update current status:

`#git pull https://github.com/MannLabs/single-cell-DVP.git`

Push changes by indicating the date of change and editor (e.g.: Florian Rosenberger on Feb 20 e.g. 20230220_FR)

1. Get the path of you local GitHub folder

`P:\03_Experiments\24_Borderline_Project\18_Github\Borderline_Manuscript>git init`

2. Add files and commit changes:

`>git add .`

`>git commit -m "20230220_Editor"`

`>git remote -v`

`>git push origin main`



