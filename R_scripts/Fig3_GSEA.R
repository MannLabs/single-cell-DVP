###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## Read relevant data
load("../output/variables/meta_distances_bins.R")
load("../output/Variables/d.R")
load("../output/Variables/meta_pg.R")

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
                                                   Webgestalt_group_8), type='full')
  
  save(webgestalt_combined_score, file = "../data/webgestalt_combined_score.KEGG.R")
  save(webgestalt_combined_score, file = "../output/Variables/webgestalt_combined_score.KEGG.R")
}

if(decision_webgestalt == 2){
  load("../data/webgestalt_combined_score.KEGG.R")
}

geneset_incl <- c("mmu00982", "mmu00020", "mmu00220", "mmu00190", "mmu05225", "mmu04932", "mmu03320", "mmu00140", "mmu01040")

webgestalt_combined_score %>%
  filter(geneSet %in% geneset_incl) %>%
  mutate(significant = FDR < 0.1) %>%
  ggplot(aes(x = -(cluster_peak_ID-9), y = normalizedEnrichmentScore, group = description))+
  geom_line(aes(color=description), size = 2)+
  geom_point(aes(size = significant), alpha = 0.3)+
  scale_color_manual(values = viridis(11))+
  facet_wrap(.~description)+
  theme_classic()+
  geom_hline(yintercept = 0, lty = "dotted") -> plot_pathwayEnrichment

ggsave(plot_pathwayEnrichment, file = "../output/Figures/pathwayEnrichment.pdf" , width = 10, height = 7)
