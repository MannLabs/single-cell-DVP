###########################
#### scDVP Figure Code ####
###########################

#### -- Data wrangling -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## -- dditional function
source("./Functions/normalize_core.R")
load("../output/Variables/refquant.R")

## -- Read common contaminants
contaminants <- read_tsv("../data/external/221109_common-contaminants_cRAP.tsv")

## -- Read image to MS crossreference file
img_meta <- read_csv("../data/meta/meta_img-proteome.csv") %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))

## -- Read fluorescent instensity values per shape
img_fluovalues <- read_csv("../data/meta/meta_img-quantification.csv") %>%
  dplyr::rename(bio_ID = Slide) %>%
  left_join(img_meta) %>%
  drop_na(cell_ID) %>%
  mutate(Area = Area*0.0841/100) # 0.29^2 = 0.0841 is the scaling factor from pxl to um on the microscope, divided by 100 which is a default in BIAS

## -- Read geometric distance data
meta_distances <- read_csv("../data/meta/meta_img-geometric-distance.csv") %>%
  mutate(sum_distance = `dist_to_portal_vein(pixels)` + `dist_to_central_vein(pixels)`) %>%
  mutate(ratio = `dist_to_portal_vein(pixels)` / sum_distance) %>%
  drop_na(ratio) %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))

write_tsv(meta_distances, file = "../output/Tables/scDVP-geometric-distances.tsv")

## -- Read proteomics data
label_pattern = "\\(Dimethyl-n-[0-9]+\\)"

read_tsv(refquant) %>%
  gather(R.FileName, int, !protein) %>%
  # replace(is.na(.), 0) %>%
  # drop_na(int) %>%
  #mutate(int = 2^as.numeric(int)) %>%
  mutate(run_ID = paste(str_replace_all(str_replace_all(R.FileName, ".*scDVP_", ""), "_.*", ""),
                        str_replace_all(str_replace_all(R.FileName, ".*AID8_", ""), "_.*", ""),
                        sep = "_")) %>%
  mutate(label = str_replace(R.FileName, ".*_", "")) %>%
  mutate(label = str_replace_all(label, "\\(", ""), label = str_replace_all(label, "\\)", "")) %>%
  mutate(cell_ID = paste(run_ID, label, sep = "_")) %>%
  mutate(bio_ID = str_replace(cell_ID, "_.*", "")) %>%
  mutate(tube_ID = str_replace(run_ID, ".*_", "")) %>%
  dplyr::rename(Protein = protein) -> d_full_set

# if(median(d_full_set$int, na.rm = T) < 30){
#   d_full_set %>%
#     mutate(int = 2^int) -> d_full_set
# }

d_full_set %>%
  filter(bio_ID %in% c("m3B", "m4A", "m5C")) %>%
  filter(!Protein %in% contaminants$Entry) %>%
  left_join(img_meta) %>%
  mutate(heps = Index <= 302) %>%
  filter(!cell_ID %in% c("m3B_14_target8", "m3B_40_target4", "m4A_58_target8", "m4A_67_target4")) -> d_all

meta_pg <- data.frame(Protein = unique(d_all$Protein)) %>%
  mutate(Symbol = mapIds(org.Mm.eg.db,
                         keys=str_replace_all(Protein, ".*;", ""),
                         column="SYMBOL",
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  mutate(Genename = mapIds(org.Mm.eg.db,
                           keys=str_replace_all(Protein, ".*;", ""),
                           column="GENENAME",
                           keytype="UNIPROT",
                           multiVals="first")) %>%
  mutate(ENSEMBL = mapIds(org.Mm.eg.db,
                          keys=str_replace_all(Protein, ".*;", ""),
                          column="ENSEMBL",
                          keytype="UNIPROT",
                          multiVals="first"))

## --  Data filtering
statTable_all <- d_all %>%
  left_join(meta_pg) %>%
  drop_na(int)

statTable_all %>%
  filter(label != "Dimethyl-n-0") %>%
  group_by(run_ID, label, cell_ID, bio_ID) %>%
  summarise(n = n()) %>%
  left_join(img_fluovalues %>% filter(Channel == "CFP")) -> statTable_all_n

CIn_upper_all <- median(statTable_all_n$n) + 3* sd(statTable_all_n$n)
CIn_lower_all <- median(statTable_all_n$n) - 1.5* sd(statTable_all_n$n)
CIArea_upper_all <- median(statTable_all_n$Area, na.rm = T) + 3* sd(statTable_all_n$Area, na.rm = T)

statTable_all_n %>%
  filter(n > CIn_lower_all, n < CIn_upper_all, (Area < CIArea_upper_all | is.na(Area))) %>%
  pull(cell_ID) -> SA_incl_all

d_all %>% 
  filter(label != "Dimethyl-n-0") %>%
  filter(cell_ID %in% SA_incl_all) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_all

## Normalization to core proteome
d_all_norm <- normalize_core(data = d_all, proteome_core = proteome_complete_all, SA_incl = SA_incl_all) %>%
  left_join(meta_pg)

## Normalizing hepatocyte data
d_heps <- d_all %>%
  filter(heps == TRUE)

d_heps %>% 
  filter(label != "Dimethyl-n-0") %>%
  filter(cell_ID %in% SA_incl_all) %>%
  dplyr::select(Protein, int, cell_ID) %>%
  spread(cell_ID, int) %>%
  filter(complete.cases(.)) %>%
  pull(Protein) -> proteome_complete_heps

d <- normalize_core(data = d_heps, proteome_core = proteome_complete_heps, SA_incl = SA_incl_all) %>%
  left_join(meta_pg)

SA_incl_heps <- unique(d$cell_ID)

## save variables to file
save(d_all_norm, file = "../Output/Variables/d_all_norm.R")
save(statTable_all, file = "../Output/Variables/statTable_all.R")
save(SA_incl_all, file = "../Output/Variables/SA_incl_all.R")
save(CIn_lower_all, file = "../Output/Variables/CIn_lower_all.R")
save(meta_pg, file = "../Output/Variables//meta_pg.R")
save(statTable_all_n, file = "../Output/Variables/statTable_all_n.R")
save(d_all, file = "../Output/Variables/d_all.R")
save(meta_distances, file = "../Output/Variables/meta_distances.R")
save(img_meta, file = "../Output/Variables/img_meta.R")
save(img_fluovalues, file = "../Output/Variables/img_fluovalues.R")
save(d_full_set, file = "../Output/Variables/d_full_set.R")
save(d, file = "../Output/Variables/d.R")
