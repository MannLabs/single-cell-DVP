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

- [Data - Supplementary_Figure_6](/data/Supplementary_Figure_6)

- (Substructures can be created in this folder. Once this is done, I can create a more detailed overview here for you)

## Results

Figures and result dataframes are saved in the [output](/output/) folder. 

- [Results - Supplementary_Figure_6](/output/Figures/Supplementary_Figure_6)

- (Substructures can be created in this folder. Once this is done, I can create a more detailed overview here for you)

## R scripts

### Data wrangling scripts
- [Filtering criteria](R_scripts/Data-wrangling.R)
- [Prediction proteome calculation](R_scripts/Prediction_class_proteomes.R)

### Scripts for main figures
#### Figure 2, Depth of single shape proteomes and estimation of nuclear compartment
- [Figure 2, Intensity of top four histones](R_scripts/Histone-levels.R)

#### Figure 3, Single shape proteomes are accurate descriptors of zonated hepatocytes
- [Figure 3a, PCA with distance metric overlay](R_scripts/PCA_vs_geometric-distance.R)
- [Figure 3b, Distance metric versus PC1](R_scripts/PCA_vs_geometric-distance.R)
- [Figure 3c, Global heatmap of 90% complete proteome](R_scripts/Heatmap_global_distances.R)
- [Figure 3d, Expression by location](R_scripts/Heatmap_markers.R)
- [Figure 3e, Relative expression of top-10 differential proteins](R_scripts/Spatial_expression_top10.R)
- [Figure 3f, Spatial Gene Set Enrichment Analysis]()
- [Figure 3g, Urea cycle and peroxisomal fatty acid metabolism](R_scripts/Pathway_Urea_Peroxisome.R)

#### Figure 4, Combining imaging and proteome data for a machine-learned model
- [Figure 4a, Pseudo-FACS plot](R_scripts/Pseudo-FACS.R)
- [Figure 4b, Proteome bin versus staining intensity](R_scripts/Pseudo-FACS.R)
- [Figure 4e, Proteome prediction of an unseen section](R_scripts/Prediction_new_mouse.R)

### Scripts for supplementary figures
#### Supplementary Figure S3: Performance overview of single shape proteomes
- [Figure S3a, Labelling efficiency](R_scripts/Labelling-efficiency.R)
- [Figure S3b, Protein IDs per biological replicate](R_scripts/Protein-IDs_vs_Runs.R)
- [Figure S3c, Abundance range and transcription factors](R_scripts/ank_versus_Intensity.R)
- [Figure S3d, Microdissected area versus protein IDs](R_scripts/Protein-IDs_vs_Area.R)
- [Figure S3e, Data completeness versus intensity](R_scripts/Completeness_vs_Intensity.R)
- [Figure S3f, Normalization versus CVs](R_scripts/CVs.R)

#### Supplementary Figure S4: Concatentation of shapes
- [Figure S4](R_scripts/PCA_reductive.R)

#### Supplementary Figure S5: Dimensionality reduction of single shape data
- [Figure S4a, PCA with Ass1 expression levels](R_scripts/PCA_Hepatocytes.R)
- [Figure S4b, PCA with Cyp2e1 expression levels](R_scripts/PCA_Hepatocytes.R)
- [Figure S4c, Distance metric versus PC12](R_scripts/PCA_vs_geometric-distance.R)
- [Figure S4d, PCA with loadings](R_scripts/PCA_Hepatocytes.R)
- [Figure S4e, PCA including endothelial cells](R_scripts/PCA_Endothelial.R)

#### Supplementary Figure S6: Normality checks of single shape data
- [Figure S5a, Volcano plot]()
- [Figure S5b, Shapiro-Wilk-Test](R_scripts/Shapiro.R)
- [Figure S5c, Pearson correlation between zones](R_scripts/)
- [Figure S5d, Relative expression of bottom-10 differential proteins](R_scripts/Shapiro.R)

#### Supplementart Figure S7
- [Script - Supplementary_Figure_6](R_scripts/Supplementary_Figure_6.R)

- [Template](Template.R)

- An exemplary script is saved. We could keep all headers of the R scripts consistent, thereby providing a short overview of the aim and content of each file. 
- Additonally, the header contains a function that is referring the folder structure of each script to the respective local copies of the GitHub folder. This means that anyone who will download the folder can run the scripts. Here, the data will be taken from the data folder and the output will be saved in the output folder. 

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



