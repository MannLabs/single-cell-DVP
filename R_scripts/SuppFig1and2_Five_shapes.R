###########################
#### scDVP Figure Code ####
###########################

#### -- Figure XX -- ####

## -- Prepare Workspace
cat("\014")
rm(list=ls())

## -- Read data
d <- read_tsv("../data/protein/proteintable_fiveshape.tsv") %>%
  drop_na(Protein.Names) %>%
  dplyr::select(-Genes, -Protein.Ids, -Protein.Names, -First.Protein.Description) %>%
  gather(file, value, !Protein.Group) %>%
  mutate(value = log2(value)) %>%
  mutate(sample = str_replace_all(file, ".*DIA_", "")) %>%
  mutate(sample = str_replace_all(sample, "_S.*", "")) %>%
  dplyr::select(-file)

write_tsv(d %>% spread(sample, value), file = "../output/Tables/Five-shape-proteome.tsv")

read_tsv("../data/protein/proteintable_fiveshape.tsv") %>%
  drop_na(Protein.Names) %>%
  dplyr::select(Genes, Protein.Ids, Protein.Names, First.Protein.Description) -> d_names

## -- Data filtering
# Chose sample to include based on distribution of number of protein (+/- 1.5 SD)

d %>%
  group_by(sample) %>%
  drop_na() %>%
  summarise(n = n()) -> d_n

CIn_upper <- median(d_n$n) + 1.5* sd(d_n$n)
CIn_lower <- median(d_n$n) - 1.5* sd(d_n$n)

d_n %>%
  filter(n > CIn_lower, n < CIn_upper) %>%
  pull(sample) -> SA_incl

# Subset based on =/- 1.5SD, update data and meta
meta <- read_tsv("../data/meta/meta_fiveshape.txt") %>%
  filter(sample %in% SA_incl) %>%
  arrange(sample) %>%
  column_to_rownames("sample")

read_tsv("../data/meta/meta_fiveshape.txt") %>%
  mutate(included = sample %in% SA_incl) %>%
  arrange(sample) %>%
  write_tsv(., file = "../output/Tables/Five-shape-meta.tsv")

d %>%
  filter(sample %in% SA_incl) %>%
  spread(sample, value) %>%
  column_to_rownames("Protein.Group") -> d_sub

## -- QC of samples
d %>%
  drop_na(value) %>%
  group_by(sample) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(`Sample ID` = c(1:length(unique(sample)))) %>%
  ggplot(aes(x = `Sample ID`, y = n))+
  geom_smooth(se = F)+
  geom_point(pch = 21, color = "cornflowerblue", size = 2) +
  scale_y_continuous(limits = c(0,2000))+
  theme_classic() +
  geom_hline(yintercept = CIn_upper, lty = "dotted") +
  geom_hline(yintercept = CIn_lower, lty = "dotted") -> plot

d %>%
  drop_na(value) %>%
  group_by(sample) %>%s
  summarise(n = n()) %>%
  summarise(mean = mean(n), median = median(n))

## -- PCAs

# Use complete cases for PCA
d_sub %>%
  filter(complete.cases(.)) -> d_complete

# Compute four clusters by kmeans clustering. Biological rational with common portal vein, periportal, pericentral, central vein distinction
clusters_compl <- kmeans(t(d_complete), centers = 4, iter.max = 1000, nstart = 50)

as.data.frame(clusters_compl[["cluster"]]) %>%
  rownames_to_column("ID") -> cluster_compl_df

colnames(cluster_compl_df)[2] = "cluster"

meta %>%
  rownames_to_column("ID") %>%
  left_join(cluster_compl_df) %>%
  column_to_rownames("ID") -> meta_cluster

p_cluster <- PCAtools::pca(d_complete, metadata = meta_cluster, removeVar = 0.1)

write_tsv(as.data.frame(p_cluster$rotated) %>% rownames_to_column("sample"), file = "../output/Tables/Five-shape-PC.tsv")

PCAtools::biplot(p_cluster ,
                 colby = 'cluster',
                 colLegendTitle = 'Cluster',
                 # encircle config
                 encircle = TRUE,
                 encircleFill = TRUE,
                 hline = 0, vline = c(-25, 0, 25),
                 legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0,
                 showLoadings = T, lab = NA)+
  scale_colour_viridis() -> plot_pca_kmeans

## Asl
meta %>%
  rownames_to_column("ID") %>%
  left_join(cluster_compl_df) %>%
  left_join(d %>%
              filter(Protein.Group == "Q91YI0") %>%
              dplyr::rename(ID = sample)) %>%
  column_to_rownames("ID")-> meta_cluster_asl

p_cluster <- PCAtools::pca(d_complete, metadata = meta_cluster_asl, removeVar = 0.1)

PCAtools::biplot(p_cluster ,
              colby = 'value',
              colLegendTitle = 'Cluster',
              # # encircle config
              encircle = FALSE,
              encircleLineSize = 3,
              encircleLineCol = "grey40",
              encircleFill = FALSE,
              hline = 0, vline = c(-25, 0, 25),
              legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0, lab = NULL) +
         scale_color_viridis(option = "inferno") -> plot_pca_asl

## Cyp2e1
meta %>%
  rownames_to_column("ID") %>%
  left_join(cluster_compl_df) %>%
  left_join(d %>%
              filter(Protein.Group == "Q05421") %>%
              dplyr::rename(ID = sample)) %>%
  column_to_rownames("ID")-> meta_cluster_cyp2e1

p_cluster <- PCAtools::pca(d_complete, metadata = meta_cluster_cyp2e1, removeVar = 0.1)

PCAtools::biplot(p_cluster ,
              colby = 'value',
              colLegendTitle = 'Cluster',
              # # encircle config
              encircle = FALSE,
              encircleLineSize = 3,
              encircleLineCol = "grey40",
              encircleFill = FALSE,
              hline = 0, vline = c(-25, 0, 25),
              legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0, lab = NULL) +
         scale_color_viridis(option = "inferno") -> plot_pca_cyp2e1

## --- Boxplots by marker proteins
markers <- data.frame(Genes = c("Glul", "Cyp2e1", "Ass1", "Asl", "Cyp2f2", "Cps1", "Pck1", "Cyp1a2"),
                      type = c("CV", "CV", "PV", "PV", "PV", "PV", "PV", "CV"))

d_names %>%
  filter(Genes %in% markers$Genes) %>%
  pull(Protein.Ids) -> markers_Ids

d %>%
  left_join(meta_cluster %>% rownames_to_column("sample")) %>%
  drop_na(cluster) %>%
  left_join(d_names, by = c("Protein.Group" = "Protein.Ids")) %>%
  filter(Protein.Group %in% markers_Ids) %>%
  dplyr::select(Genes, cluster, value) %>%
  left_join(markers) %>%
  group_by(Genes, cluster, type) %>%
  summarise(median = median(value, na.rm = T)) %>%
  ungroup() %>%
  group_by(Genes) %>%
  mutate(rank = rank(median)) %>%
  mutate(rank_corrected = ifelse(type == "CV", rank, 9 - rank)) %>%
  group_by(cluster) %>%
  summarise(true_cluster_mean = mean(rank_corrected)) %>%
  mutate(true_cluster = rank(true_cluster_mean)) -> true_clusters

d %>%
  left_join(meta_cluster %>% rownames_to_column("sample")) %>%
  drop_na(cluster) %>%
  left_join(d_names, by = c("Protein.Group" = "Protein.Ids")) %>%
  filter(Protein.Group %in% markers_Ids) %>%
  dplyr::select(Genes, cluster, value) %>%
  left_join(markers) %>%
  left_join(true_clusters) %>%
  dplyr::select(Genes, value, true_cluster) %>%
  group_by(Genes, true_cluster) %>%
  summarise(median = median(value, na.rm = T)) %>%
  arrange(true_cluster) %>%
  spread(Genes, median) %>%
  column_to_rownames("true_cluster") -> marker_heat_in

marker_heat_in_z <- scale(marker_heat_in)
pheatmap(t(marker_heat_in_z), cluster_rows = T, cluster_cols = F, color = inferno(50), cellwidth = 20, cellheight = 20) -> plot_heat_marker

## -- ANOVA, Volcano and Heatmap

# Sample filtering
d_sub %>%
  filter(rowSums(is.na(.)) < 0.8 * length(colnames(d_sub))) %>%
  rownames_to_column("Protein_ID") %>%
  gather(sample, value, !Protein_ID) %>%
  filter(!sample %in% c("slide01_frame20", "slide10_frame12", "slide13_frame06")) %>%
  spread(sample, value) %>%
  column_to_rownames("Protein_ID") -> d_sub_outliers

# Clusters have to be ordered first based on marker genes defined above
meta_cluster_true <- meta_cluster %>%
  rownames_to_column("proteome_name") %>%
  left_join(true_clusters) %>%
  filter(proteome_name %in% colnames(d_sub_outliers))

# Stats into ANOVA
design <- cbind(Grp1 = 1, as.factor(meta_cluster_true[match(colnames(d_sub_outliers), meta_cluster_true$proteome_name),]$true_cluster))

fit <- lmFit(d_sub_outliers, design)
fit <- eBayes(fit)
d_anova_out <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  mutate(Protein.Ids = rownames(.)) %>%
  left_join(d_names) %>%
  mutate(EnsemblGeneID = mapIds(org.Mm.eg.db, 
                                keys=Protein.Ids, 
                                column="ENSEMBL", 
                                keytype="UNIPROT",
                                multiVals="first"))

#write_csv(d_anova_out, file = "./Supp1_FiveShape_ANOVA.tsv")

length(d_anova_out %>%
         filter(adj.P.Val < 0.01) %>% pull (Protein.Ids))

# Volcano plot
ggplot() +
         geom_point(data = d_anova_out, aes(y = -log10(adj.P.Val), x = logFC), alpha = 0.5, size = 1.5) +
         geom_point(data = d_anova_out %>% filter(adj.P.Val < 0.01 & logFC > 0.5), aes(y = -log10(adj.P.Val), x = logFC), alpha = 0.8, color = "darkorange", size = 2)+
         geom_point(data = d_anova_out %>% filter(adj.P.Val < 0.01 & logFC < -0.5), aes(y = -log10(adj.P.Val), x = logFC), alpha = 0.8, color = "darkblue", size = 2)+
         geom_point(data = d_anova_out %>% filter(adj.P.Val < 0.01 & logFC < 0.5 & logFC > 0), aes(y = -log10(adj.P.Val), x = logFC), alpha = 0.4, color = "darkorange", size = 2)+
         geom_point(data = d_anova_out %>% filter(adj.P.Val < 0.01 & logFC > -0.5 & logFC < 0), aes(y = -log10(adj.P.Val), x = logFC), alpha = 0.4, color = "darkblue", size = 2)+
         theme_bw()+
         geom_hline(yintercept = -log10(0.01), lty = "dotted")+
         geom_text_repel(data = d_anova_out %>% filter(adj.P.Val < 0.01 & abs(logFC) > 0.5), aes(y = -log10(adj.P.Val), x = logFC, label = Genes))+
         geom_vline(xintercept = 0.5, lty = "dotted") +
         geom_vline(xintercept = -0.5, lty = "dotted") -> plot_volcano

d_anova_out %>%
  dplyr::select(Protein.Ids, logFC) %>%
  arrange(logFC) -> gsea_in

## -- GSEA
d_anova_out %>%
  filter(adj.P.Val < 0.05) %>%
  filter(logFC > 0) %>%
  pull(Protein.Ids) -> proteins_clust_cv

d_anova_out %>%
  filter(adj.P.Val < 0.05) %>%
  filter(logFC < 0) %>%
  pull(Protein.Ids) -> proteins_clust_pv

WebGestaltR_ora_cv <- WebGestaltR(enrichMethod="ORA", organism="mmusculus",
                                  enrichDatabase = c("pathway_Wikipathway"),
                                  interestGene = proteins_clust_cv,
                                  interestGeneType ="uniprotswissprot",
                                  referenceGene = d_anova_out$Protein.Ids,
                                  referenceGeneType="uniprotswissprot") %>%
  arrange(-log10(FDR)) %>%
  mutate(description = factor(description, levels = description))

WebGestaltR_ora_pv <- WebGestaltR(enrichMethod="ORA", organism="mmusculus",
                                  enrichDatabase = c("pathway_Wikipathway"),
                                  interestGene = proteins_clust_pv,
                                  interestGeneType ="uniprotswissprot",
                                  referenceGene = d_anova_out$Protein.Ids,
                                  referenceGeneType="uniprotswissprot")%>%
  arrange(-log10(FDR)) %>%
  mutate(description = factor(description, levels = description))


ggplot(data = WebGestaltR_ora_cv, aes(y = description, x = -log10(FDR + 0.00001)))+
  geom_bar(stat="identity", fill = "grey20")+
  theme_classic() -> plot_ora_cv

ggplot(data = WebGestaltR_ora_pv, aes(y = description, x = -log10(FDR + 0.00001)))+
  geom_bar(stat="identity", fill = "grey80", color = "black")+
  theme_classic() -> plot_ora_pv

# Heatmap
d_sub_outliers %>%
  rownames_to_column("Protein.Group") %>%
  filter(Protein.Group %in% (d_anova_out %>% filter(adj.P.Val <= 0.05) %>% pull(Protein.Ids))) %>%
  column_to_rownames("Protein.Group") -> d_significant

as.data.frame(p_cluster$rotate) %>%
  rownames_to_column("slide") %>%
  filter(slide %in% colnames(d_significant)) %>%
  dplyr::select(slide, PC1) %>%
  arrange(PC1) %>%
  pull(slide) -> rank_slide

d_anova_out %>%
  filter(adj.P.Val <= 0.05) %>%
  arrange(logFC) %>%
  pull(Protein.Ids) -> rank_protein

mx_significant <- scale(t(as.matrix(d_significant)))

mx_significant <- mx_significant[rank_slide, rank_protein]

myBreaks <- c(seq(-2.5,2.5, by = 0.1))
myColor <- colorRampPalette(c("#0C80C0", "#FFFFFF", "#ED2224"))(length(myBreaks))

pheatmap(mx_significant, breaks = myBreaks, color = myColor, show_colnames = T, show_rownames = F, cluster_cols = F, cluster_rows = F,
         cellheight = 1, cellwidth = 1, gaps_col = length(proteins_clust_pv)) -> p_heatmap

## -- Plot figures
ggsave(plot_pca_cyp2e1, file = "../output/Figures/Five-shape_PCA_Cyp2e1.pdf", width = 5, height = 5)
ggsave(plot_pca_asl, file = "../output/Figures/Five-shape_PCA_Asl.pdf", width = 5, height = 5)
ggsave(plot_pca_kmeans, file = "../output/Figures/Five-shape_PCA_kmeans.pdf", width = 5, height = 5)
ggsave(plot_volcano, file = "../output/Figures/Five-shape_Volcano.pdf", width = 8, height = 8)
ggsave(plot_ora_cv, file = "../output/Figures/Five-shape_ORA_CV.pdf", width = 4, height = 3)
ggsave(plot_ora_pv, file = "../output/Figures/Five-shape_ORA_PV.pdf", width = 4, height = 3)
ggsave(p_heatmap, file = "../output/Figures/Five-shape_heatmap.pdf", width = 8, height = 8)

## -- Write table
write_tsv(d_anova_out, file = "../output/Tables/Five-shape-ANOVA.tsv")
write_tsv(WebGestaltR_ora_cv, file = "../output/Tables/Five-shape-ORA-CV.tsv")
write_tsv(WebGestaltR_ora_pv, file = "../output/Tables/Five-shape-ORA-PV.tsv")
