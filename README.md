Florian Rosenberger et al., 2023 - in revision

# Spatial single-cell mass spectrometry defines zonation of the hepatocyte proteome

Here could be a nice staining picture or a graphical abstract

## Abstract

Single-cell proteomics by mass spectrometry (MS) is emerging as a powerful and unbiased method for the characterization of biological heterogeneity. So far, it has been limited to cultured cells, whereas an expansion of the method to complex tissues would greatly enhance biological insights. Here we describe single-cell Deep Visual Proteomics (scDVP), a technology that integrates high-content imaging, laser microdissection and multiplexed MS. scDVP resolves the context-dependent, spatial proteome of murine hepatocytes at a current depth of 1,700 proteins from a slice of a cell. Half of the proteome was differentially regulated in a spatial manner, with protein levels changing dramatically in proximity to the central vein. We applied machine learning to proteome classes and images, which subsequently inferred the spatial proteome from imaging data alone. scDVP is applicable to healthy and diseased tissues and complements other spatial proteomics or spatial omics technologies.

## Table of contents

1. [Data repository](#Data-repository)
2. [Results](#Results)
3. [R Scripts](#R-Scripts)
4. [GitHub Notes](#GitHub-Notes)  

## Data repository

Processed mass spectrometry raw data and other input files have been saved in the following folders:

- [Data - Supplementary_Figure_6](/data/Supplementary_Figure_6)

- (Substructures can be created in this folder. Once this is done, I can create a more detailed overview here for you)

## Results

Figures and result dataframes are saved in the [output](/output/) folder. 

- [Results - Supplementary_Figure_6](/output/Supplementary_Figure_6)

- (Substructures can be created in this folder. Once this is done, I can create a more detailed overview here for you)

## R scripts

[...]

- [Script - Supplementary_Figure_6](Supplementary_Figure_6.R)

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



