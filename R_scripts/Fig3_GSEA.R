###########################
#### scDVP Figure Code ####
###########################

#### -- Figure 3F -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/meta_distances.R")
load("../output/Variables/d.R")
load("../output/Variables/meta_pg.R")

## Define number of classes
classes = 20

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> meta_distances_bins

meta_distances_bins %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bins) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> meta_distances_bins

## Filter on 70% complete proteins 
d %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  column_to_rownames("Protein") %>%
  filter(rowSums(is.na(.)) / ncol(.) <= 0.5) -> d_wide_70
  
decision_webgestalt <- 1

if(decision_webgestalt == 1){
  for(i in c(1:max(meta_distances_bins$bin))){
    
    meta_distances_bins %>%
      mutate(GOI = bin == i) -> meta_tmp
    
    d_extreme <- d_wide_70[,rownames(meta_tmp)]
    
    design_extreme <- model.matrix(~meta_tmp[colnames(d_extreme),]$GOI)
    
    fit <- lmFit(d_extreme, design_extreme)
    fit <- eBayes(fit)
    limma_anchor_goi <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
      rownames_to_column("Protein") %>%
      left_join(meta_pg)
    
    limma_anchor_goi %>%
      dplyr::select(ENSEMBL, logFC, adj.P.Val) %>%
      mutate(score = logFC * -log10(adj.P.Val)) %>%
      dplyr::select(ENSEMBL, score) -> gsea_i
    
    WebGestaltR(
      enrichMethod = "GSEA",
      organism = "mmusculus",
      enrichDatabase = c("pathway_KEGG"),
      interestGene = gsea_i,
      interestGeneType = "ensembl_gene_id",
      fdrThr = 1.0,
      isOutput = F) %>%
      mutate(cluster_peak_ID = i) -> webgestalt_out
    
    goi_webgestalt_name <- paste("Webgestalt_group", i, sep = "_")
    goi_limma_name <- paste("limma_group", i, sep = "")
    
    assign(goi_limma_name, limma_anchor_goi)
    assign(goi_webgestalt_name, webgestalt_out)
  }
  
  webgestalt_combined_score <- plyr::join_all(list(Webgestalt_group_1,
                                                   Webgestalt_group_2,
                                                   Webgestalt_group_3,
                                                   Webgestalt_group_4,
                                                   Webgestalt_group_5,
                                                   Webgestalt_group_6,
                                                   Webgestalt_group_7,
                                                   Webgestalt_group_8,
                                                   Webgestalt_group_9,
                                                   Webgestalt_group_10,
                                                   Webgestalt_group_11,
                                                   Webgestalt_group_12,
                                                   Webgestalt_group_13,
                                                   Webgestalt_group_14,
                                                   Webgestalt_group_15,
                                                   Webgestalt_group_16,
                                                   Webgestalt_group_17,
                                                   Webgestalt_group_18,
                                                   Webgestalt_group_19,
                                                   Webgestalt_group_20), type='full')
  
  save(webgestalt_combined_score, file = "../data/webgestalt_combined_score.KEGG.R")
  save(webgestalt_combined_score, file = "../output/Variables/webgestalt_combined_score.KEGG.R")
}

  if(decision_webgestalt == 2){
  load("../data/webgestalt_combined_score.KEGG.R")
}

geneset_incl <- c("mmu00982", "mmu00220", "mmu00190", "mmu03320", "mmu00140", "mmu01212")

webgestalt_combined_score %>%
  mutate(FDR_manual = p.adjust(pValue, method = "fdr")) %>%
  filter(geneSet %in% geneset_incl) %>%
  mutate(significant = FDR_manual < 0.1) %>%
  ggplot(aes(x = -(cluster_peak_ID-classes-1), y = enrichmentScore, group = description))+
  geom_line(aes(color=description), size = 2)+
  geom_point(aes(size = significant), alpha = 0.3)+
  scale_color_manual(values = viridis(6)[c(2,3,6,4,5,1)])+
  facet_wrap(.~description)+
  theme_classic()+
  geom_hline(yintercept = 0, lty = "dotted") -> plot_pathwayEnrichment

webgestalt_combined_score %>%
  mutate(FDR_manual = p.adjust(pValue, method = "fdr")) %>%
  filter(FDR_manual < 0.1) %>%
  distinct(geneSet) %>%
  pull(geneSet) -> geneSet_significant

webgestalt_combined_score %>%
  mutate(FDR_manual = p.adjust(pValue, method = "fdr")) %>%
  filter(geneSet %in% geneSet_significant) %>%
  group_by(geneSet) %>%
  arrange(cluster_peak_ID) %>%
  mutate(diff = normalizedEnrichmentScore - lead(normalizedEnrichmentScore)) %>%
  slice_max(abs(diff)) %>%
  mutate(bin = -(cluster_peak_ID-classes-1)) %>%
  dplyr::select(geneSet, description, bin, diff) -> test

ggsave(plot_pathwayEnrichment, file = "../output/Figures/pathwayEnrichment.pdf" , width = 11, height = 7)

## -- Write tables
write_tsv(webgestalt_combined_score, file = "../output/Tables/scDVP_GSEA.tsv")
