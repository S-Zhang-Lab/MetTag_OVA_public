### Step7.3.3 analysis UTSW02 - immune cell focused ###

# Load Packages + DATA ###
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(scales)
library(plotly)

# Set seed
set.seed(42)

#set working directory for loading in raw data
setwd("/Users/s438978/Desktop/RData/MetTag_U6_LARRY_BC_OVA")
setwd("C:/Users/Zhang Lab/Desktop/Emma/R01_and_2_w_LARRY")

#Load in data
#Note: later SCT versions with date at the end will work - transcriptome was not changed
# just LARRY BC assignment 
Ova_merged_SCT <- readRDS("Ova_merged_SCT_EA.rds")

# 3_ IMMUNE RELATED ANALYSIS
### Go back to large object containing all cells and look into innate immune 
DimPlot(Ova_merged_SCT, label = TRUE)
Ova_merged_SCT@meta.data$Broad_cell.ID <- recode(Ova_merged_SCT@meta.data$seurat_clusters, "36" = "Mast", "42" = "pDC", "0" = "ID8.Tumor", "1" = "ID8.Tumor", "4" = "ID8.Tumor", "10" = "ID8.Tumor", "11" = "ID8.Tumor", "24" = "ID8.Tumor", "26" = "ID8.Tumor", "38" = "ID8.Tumor", "41" = "ID8.Tumor", "5" = "Myeloid", "7" = "Myeloid", "8" = "Myeloid", "16" = "Myeloid", "18" = "Myeloid", "19" = "Myeloid", "29" = "Myeloid", "27" = "Myeloid", "32" = "Myeloid", "34" = "Myeloid", "33" = "Myeloid", "37" = "Myeloid", "2" = "B.cell", "3" = "B.cell", "9" = "B.cell", "14" = "B.cell", "12" = "B.cell", "22" = "B.cell", "28" = "B.cell", "40" = "B.cell", "42" = "B.cell", "13" = "T/NK.cell", "15" = "T/NK.cell", "17" = "T/NK.cell", "21" = "T/NK.cell", "23" = "T/NK.cell", "30" = "T/NK.cell", "31" = "T/NK.cell", "6" = "Stromal", "20" = "Stromal", "25" = "Stromal", "35" = "Stromal", "39" = "Stromal")
Idents(Ova_merged_SCT) <- "Broad_cell.ID"
DimPlot(Ova_merged_SCT, label = TRUE)

#Compare Cluster Proportions Between Control and Experimental 
Ova_merged_SCT@meta.data$Combined.HTO_group.order <- recode(Ova_merged_SCT@meta.data$Combined.HTO_group, "asc:3_Lavage" = "01_Lavage", "asc:1_AscMet" = "02_AscMet", "om:4_NaiveOm" = "03_NaiveOm", "om:2_OmMet" = "04_OmMet")
table(Ova_merged_SCT@meta.data$Broad_cell.ID,Ova_merged_SCT@meta.data$Combined.HTO_group.order)
freq_table <- prop.table(x = table(Ova_merged_SCT@meta.data$Broad_cell.ID,Ova_merged_SCT@meta.data$Combined.HTO_group.order),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(Ova_merged_SCT@meta.data$Broad_cell.ID) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)

Idents(Ova_merged_SCT) <- "Broad_cell.ID"
Myeloids <- subset(Ova_merged_SCT, idents = c("Myeloid"))
T.NK <- subset(Ova_merged_SCT, idents = c("T/NK.cell"))
#re-cluster subsetted cell type
Idents(T.NK) <- "seurat_clusters"
T.NK <- RunPCA(T.NK, verbose = FALSE)
T.NK <- RunUMAP(T.NK, dims = 1:10, verbose = FALSE)
T.NK <- FindNeighbors(T.NK, dims = 1:30, verbose = FALSE)
T.NK <- FindClusters(T.NK, verbose = FALSE, resolution = 0.2)
DimPlot(T.NK, label = TRUE)
DimPlot(T.NK, split.by = "Combined.HTO_group.order")
#export dimplots    

#re-cluster subsetted cell type
Idents(Myeloids) <- "seurat_clusters"
Myeloids <- RunPCA(Myeloids, verbose = FALSE)
Myeloids <- RunUMAP(Myeloids, dims = 1:10, verbose = FALSE)
Myeloids <- FindNeighbors(Myeloids, dims = 1:30, verbose = FALSE)
Myeloids <- FindClusters(Myeloids, verbose = FALSE, resolution = 0.2)
DimPlot(Myeloids, label = TRUE)
DimPlot(Myeloids, split.by = "Combined.HTO_group.order")
#export dimplots 

# Check identity of increasing populations (C0,1 and C5,6)
VlnPlot(Myeloids, features = c("Adgre1", "Cd68", "Csf1r", "Itgam", "Mertk", "Trem2"))

# examine DEGs between Asc and OmMet macrophages
# GATING 
# Make two sub-populations: Adgre1+Cd68+ and Adgre1- Cd68-
Idents(Myeloids) <- "seurat_clusters"
plotAdgre1xCd68 <- FeatureScatter(Myeloids, feature1 = "Adgre1", feature2 = "Cd68")
#select double negative
Myeloids <- CellSelector(plot = plotAdgre1xCd68, object = Myeloids, ident = "Non-MAC")
#select double positive
Myeloids <- CellSelector(plot = plotAdgre1xCd68, object = Myeloids, ident = "MAC")
DimPlot(Myeloids, split.by = "Combined.HTO_group.order")
# export 4 x 4 pdf of gate

Myeloid_subset <- subset(Myeloids, idents = c("Non-MAC", "MAC"))

Macs <- subset(Myeloid_subset, ident = "MAC")
Macs <- PrepSCTFindMarkers(Macs)
Idents(Macs) <- "Combined.HTO_group"

DEG_Macs <- FindMarkers(Macs, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", verbose = TRUE)
write.csv(DEG_Macs, file = "DEG_Macs_All_Om_vs_Asc.csv")

# Load ggpubr if not already
library(ggpubr)

# Define signatures
m1_genes <- c(
  # Canonical cytokines & pro-inflammatory mediators
  "Nos2", "Tnf", "Il1b", "Il6", "Il12b", "Il18", "Cxcl9", "Cxcl10", "Cxcl11", "Ccl5", "Ccl2",
  
  # Costimulatory molecules / activation markers
  "Cd86", "Cd80", "Cd40", "Tnfrsf1b", "Tnfrsf9", "Cd83",
  
  # Pattern recognition & IFN response
  "Tlr2", "Tlr4", "Nfkb1", "Irf1", "Irf5", "Stat1", "Gbp2", "Gbp5", "Ifit1", "Ifit3", "Rsad2",
  
  # Antigen presentation / MHC II
  "H2-Ab1", "H2-Aa", "H2-Eb1", "Cd74", "Ciita",
  
  # Metabolic & oxidative
  "Gch1", "Slc2a1", "Pfkm", "Sod2", "Txnip", "Aldoa",
  
  # Surface / activation markers
  "Cd38", "Gpr18", "Fpr2", "Cx3cr1", "Ly6c2", "Ly6a",
  
  # Inflammatory transcriptional regulators
  "Nfkbia", "Socs3", "Zbp1", "Cebpb", "Cebpd"
)


m2_genes <- c(
  # Canonical/functionally validated
  "Arg1", "Mrc1", "Il10", "Chil3", "Retnla", "Ccl17", "Ccl22", "Tgfb1", "Tgfbi",
  
  # Surface receptors / immune modulators
  "Cd163", "Mrc1", "Cd209a", "Clec10a", "Clec4n", "Marco", "Trem2", "Pdcd1lg2",
  
  # Transcription factors / regulators
  "Maf", "Mafb", "Pparg", "Klf4", "Stat6", "Irf4",
  
  # Metabolic / mitochondrial
  "Fxyd2", "Fabp4", "Plin2", "Acly", "Scarb1", "Slc1a3", "Slc7a2",
  
  # ECM/wound repair/fibrosis
  "Spp1", "F13a1", "Serpine1", "Timp1", "Lpl", "Lgals3", "Gpnmb", "Sema3e",
  
  # Secreted factors / immunosuppressive
  "Vegfa", "Areg", "Bmp2", "Bmp7", "Socs1", "Csf1", "Il1rn", "Anxa1",
  
  # Tumor-associated macrophage (TAM) associated M2-like genes
  "Trem2", "Apoe", "C1qa", "C1qb", "C1qc", "Lair1", "Ccl8", "Ccl13", "Fcrls", "Ms4a7"
)


# Add scores
Macs <- AddModuleScore(Macs, features = list(m1_genes), name = "M1_Score")
Macs <- AddModuleScore(Macs, features = list(m2_genes), name = "M2_Score")

# Compare ascites vs omental
Idents(Macs) <- "Combined.HTO_group"
Macs <- subset(Macs, idents = c("om:2_OmMet", "asc:1_AscMet"))

# Generate violin plot with p-value using ggpubr
library(ggpubr)
library(ggplot2)

# M1 score
p1 <- VlnPlot(Macs, features = "M1_Score1", group.by = "Combined.HTO_group", pt.size = 0) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  ggtitle("M1 Score by Group")

# M2 score
p2 <- VlnPlot(Macs, features = "M2_Score1", group.by = "Combined.HTO_group", pt.size = 0) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  ggtitle("M2 Score by Group")

# Combine with patchwork
p1 + p2

### Check if hybrid state present
# Extract metadata including module scores and grouping info
df <- Macs@meta.data %>%
  select(M1_Score1, M2_Score1, Combined.HTO_group) %>%
  rename(M1 = M1_Score1, M2 = M2_Score1, Group = Combined.HTO_group)

# Basic biaxial scatter plot
ggplot(df, aes(x = M1, y = M2, color = Group)) +
  geom_point(alpha = 0.4, size = 1) +
  labs(
    title = "M1 vs M2 Module Scores",
    x = "M1 Score",
    y = "M2 Score"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set1")

library(ggplot2)
library(dplyr)


# Scatter plot with filled points and transparent filled ellipses
ggplot(df, aes(x = M1, y = M2, color = Group, fill = Group)) +
  geom_point(alpha = 0.8, size = 1.5, shape = 21, stroke = 0.2) +
  stat_ellipse(
    geom = "polygon", 
    alpha = 0.15,     # Transparency of fill
    color = NA        # No ellipse border, use fill only
  ) +
  stat_ellipse(
    geom = "path",    # Adds solid border around the ellipse
    size = 0.8
  ) +
  labs(
    title = "M1 vs M2 Module Scores",
    x = "M1 Score",
    y = "M2 Score"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")


# Scatter plot with filled points and transparent filled ellipses

ggplot(df, aes(x = M1, y = M2, color = Group, fill = Group)) +
  geom_point(alpha = 0.8, size = 1.5, shape = 21, stroke = 0.2) +
  stat_ellipse(
    geom = "polygon", 
    alpha = 0.15,     # Transparency of fill
    color = NA        # No ellipse border, use fill only
  ) +
  stat_ellipse(
    geom = "path",    # Adds solid border around the ellipse
    size = 0.8
  ) +
  labs(
    title = "M1 vs M2 Module Scores",
    x = "M1 Score",
    y = "M2 Score"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

library(dplyr)

# Define hybrid score (e.g., mean of M1 and M2)
df <- df %>%
  mutate(Hybrid_Score = (M1 + M2) / 2)

# Check groups
table(df$Group)

# Perform Wilcoxon test comparing Hybrid_Score between the two groups
wilcox_test <- wilcox.test(Hybrid_Score ~ Group, data = df)

# Show the p-value
wilcox_test$p.value
#[1] 5.532458e-60

VlnPlot(Macs, features = c("Axl", "Tyrobp", "Cd47", "Lgals3", "Hmox1", "Cd44", "Anxa1", "Trem2", "Spp1"))

#oxphos genes
VlnPlot(Macs, features = c("Ndufb8", "Cox5a", "Atp5f1a", "Uqcrc2"))

#stress genes
VlnPlot(Macs, features = c("Hmox1", "Anxa1", "Ddit3"))

### Check if macrophages are infiltrating or proliferating
Idents(Myeloids) <- "seurat_clusters"
VlnPlot(Myeloids, features = c("Mki67", "Top2a", "Ccr2", "Sell", "Timd4", "Gata6"))

# examine DEGs between TREM2 high and TREM2 low producing macs
# GATING 
# (1) Subset out macs from other myeloids
Idents(Myeloids) <- "seurat_clusters"
plotCd11bxF480 <- FeatureScatter(Myeloids, feature1 = "Itgam", feature2 = "Adgre1")
#select CD11b+F480+ only
Myeloids <- CellSelector(plot = plotCd11bxF480, object = Myeloids, ident = "F480+")
Macs <- subset(Myeloids, ident = "F480+")
#select Trem2 high and Trem2 low
Idents(Macs) <- "seurat_clusters"
plotTrem2xPDL1 <- FeatureScatter(Macs, feature1 = "Cd274", feature2 = "Trem2")
# first select double negative population 
Macs <- CellSelector(plot = plotTrem2xPDL1, object = Macs, ident = "Trem2_Pdl1_DN")
# next, select double positive population
Macs <- CellSelector(plot = plotTrem2xPDL1, object = Macs, ident = "Trem2_Pdl1_DP")
# export 4 x 4 pdf of gate
Macs <- PrepSCTFindMarkers(Macs)
DEG_Trem2 <- FindMarkers(Macs, ident.1 = "Trem2_Pdl1_DP", ident.2 = "Trem2_Pdl1_DN", verbose = TRUE, logfc.threshold = 0.05)

Macs_subset <- subset(Macs, idents = c("Trem2_Pdl1_DP", "Trem2_Pdl1_DN"))

Trem2_macs <- subset(Macs, ident = "Trem2_Pdl1_DP")
Trem2_macs <- PrepSCTFindMarkers(Trem2_macs)
Idents(Trem2_macs) <- "Combined.HTO_group"
VlnPlot(Trem2_macs, features = c("Mmp13", "Mmp14", "Fos", "Coq10b", "H2-T23", "Ccr5"))
DEG_Trem2_macs <- FindMarkers(Trem2_macs, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", verbose = TRUE)

# Volcano plot of DEGs
DEG_Macs$neg_adj.pVal <- -1*log(DEG_Macs$p_val_adj, 10)
pThresh <- 2
plot(DEG_Macs$avg_log2FC, DEG_Macs$neg_adj.pVal, col=ifelse(DEG_Macs$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 180), pch=ifelse(DEG_Macs$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-5, 10, by=2.5), col="lightgray")
abline(h=seq(0,180, by=50), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Macs$avg_log2FC, DEG_Macs$neg_adj.pVal, col=ifelse(DEG_Macs$neg_adj.pVal > pThresh, ifelse(DEG_Macs$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Macs$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Macs$neg_adj.pVal > pThresh, 1.4, 1.2))
# export volcano 7x6

ggplotly(ggplot(data = DEG_Macs,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Macs), colour=ifelse(DEG_Macs$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

# Volcano plot of DEGs Trem2 high vs Trem2 low
DEG_Trem2$neg_adj.pVal <- -1*log(DEG_Trem2$p_val_adj, 10)
pThresh <- 2

plot(DEG_Trem2$avg_log2FC, DEG_Trem2$neg_adj.pVal, col=ifelse(DEG_Trem2$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 50), pch=ifelse(DEG_Trem2$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-5, 5, by=2.5), col="lightgray")
abline(h=seq(0,50, by=10), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Trem2$avg_log2FC, DEG_Trem2$neg_adj.pVal, col=ifelse(DEG_Trem2$neg_adj.pVal > pThresh, ifelse(DEG_Trem2$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Trem2$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Trem2$neg_adj.pVal > pThresh, 1.4, 1.2))
# export volcano 6x5

ggplotly(ggplot(data = DEG_Trem2,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Trem2), colour=ifelse(DEG_Trem2$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

# Volcano plot of DEGs Trem2 high OmMet vs AscMet
DEG_Trem2_macs$neg_adj.pVal <- -1*log(DEG_Trem2_macs$p_val_adj, 10)
pThresh <- 2

plot(DEG_Trem2_macs$avg_log2FC, DEG_Trem2_macs$neg_adj.pVal, col=ifelse(DEG_Trem2_macs$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 5), pch=ifelse(DEG_Trem2_macs$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-5, 10, by=5), col="lightgray")
abline(h=seq(0,5, by=1), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Trem2_macs$avg_log2FC, DEG_Trem2_macs$neg_adj.pVal, col=ifelse(DEG_Trem2_macs$neg_adj.pVal > pThresh, ifelse(DEG_Trem2_macs$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Trem2_macs$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Trem2_macs$neg_adj.pVal > pThresh, 1.4, 1.2))
# export volcano 6x5

ggplotly(ggplot(data = DEG_Trem2_macs,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Trem2_macs), colour=ifelse(DEG_Trem2_macs$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

#Compare Cluster Proportions Between Control and Experimental
table(T.NK@meta.data$seurat_clusters,T.NK@meta.data$Combined.HTO_group.order)
freq_table <- prop.table(x = table(T.NK@meta.data$seurat_clusters,T.NK@meta.data$Combined.HTO_group.order),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(T.NK@meta.data$seurat_clusters) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)

table(Myeloids@meta.data$seurat_clusters,Myeloids@meta.data$Combined.HTO_group.order)
freq_table <- prop.table(x = table(Myeloids@meta.data$seurat_clusters,Myeloids@meta.data$Combined.HTO_group.order),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(Myeloids@meta.data$seurat_clusters) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)

Idents(T.NK) <- "seurat_clusters"
T.NK.markers <- FindAllMarkers(T.NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, recorrect_umi = FALSE)
T.NK.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
View(T.NK.markers)
top10_T.NK_DEG <- T.NK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(T.NK, features = top10_T.NK_DEG$gene) +scale_fill_gradientn(colors = c("blue", "white", "red")) 
write.csv(T.NK.markers,"T_NK_marker_genes_subsetted.csv")
DefaultAssay(T.NK) <- "RNA"
T.NK <- NormalizeData(T.NK)  # if not already normalized
T.NK <- FindVariableFeatures(T.NK)

DEG_TNK <- FindMarkers(
  T.NK,
  ident.1 = "04_OmMet",
  ident.2 = "02_AscMet",
  verbose = TRUE
)

write.csv(DEG_TNK, file = "DEG_TNK_Om_vs_Asc.csv")

# Volcano plot of DEGs
DEG_TNK$neg_adj.pVal <- -1*log(DEG_TNK$p_val_adj, 10)
pThresh <- 10
plot(DEG_TNK$avg_log2FC, DEG_TNK$neg_adj.pVal, col=ifelse(DEG_TNK$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 90), pch=ifelse(DEG_TNK$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-7.5, 7.5, by=2.5), col="lightgray")
abline(h=seq(0,90, by=20), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_TNK$avg_log2FC, DEG_TNK$neg_adj.pVal, col=ifelse(DEG_TNK$neg_adj.pVal > pThresh, ifelse(DEG_TNK$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_TNK$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_TNK$neg_adj.pVal > pThresh, 1.4, 1.2))
# export volcano 6x5

ggplotly(ggplot(data = DEG_TNK,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_TNK), colour=ifelse(DEG_TNK$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

Idents(Myeloids) <- "seurat_clusters"
Myeloids <- PrepSCTFindMarkers(Myeloids)
Myeloids.markers <- FindAllMarkers(Myeloids, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Myeloids.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
View(Myeloids.markers)
top10_Mets_DEG <- Myeloids.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Myeloids, features = top10_Mets_DEG$gene) +scale_fill_gradientn(colors = c("blue", "white", "red")) 
write.csv(Myeloids.markers,"Myeloids_marker_genes_subsetted.csv")


VlnPlot(Myeloids, features = c("Itgam", "Itgax", "Xcr1", "Trem2", "C1qa", "Il18bp", "Cx3cr1", "Ly6c2", "Adgre1"))

# Take a look at both met vs naive comparisons
Idents(Myeloids) <- "Combined.HTO_group.order"
DEG_Myeloids.Asc <- FindMarkers(Myeloids, ident.1 = "02_AscMet", ident.2 = "01_Lavage", verbose = TRUE)
DEG_Myeloids.Om <- FindMarkers(Myeloids, ident.1 = "04_OmMet", ident.2 = "03_NaiveOm", verbose = TRUE)
write.csv(DEG_Myeloids.Asc, file = "DEG_Myeloids_Asc.csv")
write.csv(DEG_Myeloids.Om, file = "DEG_Myeloids_Om.csv")

# Compare individual clusters among each other
Idents(Myeloids) <- "seurat_clusters"
Myeloids <- PrepSCTFindMarkers(Myeloids)
DEG_Cluster1_vs_0 <- FindMarkers(Myeloids, ident.1 = "1", ident.2 = "0", verbose = TRUE)
DEG_Cluster5_vs_6 <- FindMarkers(Myeloids, ident.1 = "5", ident.2 = "6", verbose = TRUE)
write.csv(DEG_Cluster1_vs_0, file = "DEG_Cluster1_vs_0_myeloid.csv")
write.csv(DEG_Cluster5_vs_6, file = "DEG_Cluster5_vs_6_myeloid.csv")

library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Prepare the data
# Assuming DEG_Cluster5_vs_6 is your data frame
df <- DEG_Cluster5_vs_6
df$gene <- rownames(df)
df$neg_log10_padj <- -log10(df$p_val_adj)

# 2. Categorize genes for coloring
# Using a log2FC threshold of 0.5 and p-adj threshold of 0.05 (adjust as needed)
df <- df %>%
  mutate(diffexpressed = case_when(
    avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up",
    avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down",
    TRUE ~ "No"
  ))

# 3. Create the Plot
# We'll select the top 15 most significant genes per side
top_up <- df %>% filter(diffexpressed == "Up") %>% slice_min(p_val_adj, n = 25)
top_down <- df %>% filter(diffexpressed == "Down") %>% slice_min(p_val_adj, n = 18)
label_data <- rbind(top_up, top_down)

# 2. Generate the Plot
ggplot(df, aes(x = avg_log2FC, y = neg_log10_padj, fill = diffexpressed)) +
  # Background and Grid
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "grey90"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(color = "black")
  ) +
  
  # The points (Outline with fill)
  geom_point(aes(size = neg_log10_padj), shape = 21, color = "black", stroke = 0.3, alpha = 0.8) +
  
  # Color palette: Black (Down), Grey (NS), Red (Up)
  scale_fill_manual(values = c("Down" = "black", "No" = "grey80", "Up" = "firebrick3")) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  
  # Reference lines
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", size = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.4) +
  
  # Enhanced Gene Labeling
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = "black",
    segment.size = 0.2,
    max.overlaps = 20, # Increases the limit to show more labels
    force = 2          # Pushes labels away from each other
  ) +
  
  # Title and Axis Labels
  labs(
    title = "Asc.LPMs vs. Om.LPMs",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P-value),
    fill = "Expression"
  )

# Load required libraries
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Load DEG result (if not already in memory)
# DEG_Cluster1_vs_0 <- read.csv("DEG_Cluster1_vs_0_myeloid.csv", row.names = 1)

# Add gene names as a column
DEG_Cluster1_vs_0$gene <- rownames(DEG_Cluster1_vs_0)

# Calculate -log10 adjusted p-value
DEG_Cluster1_vs_0 <- DEG_Cluster1_vs_0 %>%
  mutate(neg_log10_pval = -log10(p_val_adj))

# Set thresholds
logfc_threshold <- 0.25
pval_threshold <- 0.05

# Identify significantly up and downregulated genes
DEG_Cluster1_vs_0 <- DEG_Cluster1_vs_0 %>%
  mutate(Significance = case_when(
    p_val_adj < pval_threshold & avg_log2FC > logfc_threshold ~ "Upregulated",
    p_val_adj < pval_threshold & avg_log2FC < -logfc_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 30 up and down genes for labeling
top_up <- DEG_Cluster1_vs_0 %>%
  filter(Significance == "Upregulated") %>%
  arrange(-avg_log2FC) %>%
  head(30)

top_down <- DEG_Cluster1_vs_0 %>%
  filter(Significance == "Downregulated") %>%
  arrange(avg_log2FC) %>%
  head(30)

top_genes <- bind_rows(top_up, top_down)

# Volcano plot
volcano <- ggplot(DEG_Cluster1_vs_0, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red3", "Downregulated" = "dodgerblue3", "Not Significant" = "grey60")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3.5, max.overlaps = 100) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Cluster 1 vs Cluster 0",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "DEG Status")

# Print plot
print(volcano)


# Examine changes in neutrophils
Idents(Myeloids) <- "seurat_clusters"
Neutrophils <- subset(Myeloids, idents = c("4"))

Idents(Neutrophils) <- "Combined.HTO_group.order"
DEG_Neutrophils.Om <- FindMarkers(Neutrophils, ident.1 = "04_OmMet", ident.2 = "03_NaiveOm", verbose = TRUE, recorrect_umi = FALSE)

# Volcano plot of DEGs
DEG_Neutrophils.Om$neg_adj.pVal <- -1*log(DEG_Neutrophils.Om$p_val_adj, 10)
pThresh <- 2
DEG_Neutrophils.Om$hits <- ifelse(DEG_Neutrophils.Om$hits > pThresh, c("hits"), c("")) 
DEG_Neutrophils.Om$label <- ifelse(DEG_Neutrophils.Om$hits == "hits", row.names(DEG_Neutrophils.Om), c("")) 
plot(DEG_Neutrophils.Om$avg_log2FC, DEG_Neutrophils.Om$neg_adj.pVal, col=ifelse(DEG_Neutrophils.Om$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 22), pch=ifelse(DEG_Neutrophils.Om$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-10, 10, by=5), col="lightgray")
abline(h=seq(0,20, by=5), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Neutrophils.Om$avg_log2FC, DEG_Neutrophils.Om$neg_adj.pVal, col=ifelse(DEG_Neutrophils.Om$neg_adj.pVal > pThresh, ifelse(DEG_Neutrophils.Om$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Neutrophils.Om$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Neutrophils.Om$neg_adj.pVal > pThresh, 1.4, 1.2))

ggplotly(ggplot(data = DEG_Neutrophils.Om,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Neutrophils.Om), colour=ifelse(DEG_Neutrophils.Om$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

# Subset the Seurat object for the groups you're interested in
T.NK_subset <- subset(T.NK, Combined.HTO_group.order %in% c("02_AscMet", "04_OmMet"))

# Specify the order of the groups for plotting
group_order <- c("02_AscMet", "04_OmMet")

# Plot the VlnPlot with the specified group order and gene expression for Ifng
Idents(T.NK_subset) <- "seurat_clusters"
VlnPlot(T.NK_subset, features = c("Ifng"), group.by = "Combined.HTO_group.order", 
        idents = group_order, 
        pt.size = 1)

T.NK <- SCTransform(T.NK, verbose = FALSE)
Idents(T.NK) <- "Combined.HTO_group.order"

DEG_T.NK <- FindMarkers(T.NK, ident.1 = "04_OmMet", ident.2 = "02_AscMet", verbose = TRUE)
write.csv(DEG_T.NK, file = "DEG_T_NK_OmMet_vs_AscMet.csv")

### gate out TREM2 only - compare high and low - no PD1
#select Trem2 high and Trem2 low
Idents(Macs) <- "seurat_clusters"
plotTrem2xCD11b <- FeatureScatter(Macs, feature1 = "Itgam", feature2 = "Trem2")
# first select double negative population 
Macs <- CellSelector(plot = plotTrem2xCD11b, object = Macs, ident = "Trem2-low")
# next, select double positive population
Macs <- CellSelector(plot = plotTrem2xCD11b, object = Macs, ident = "Trem2-high")
# export 4 x 4 pdf of gate


# Prep data for DE
Macs <- PrepSCTFindMarkers(Macs)

# Run differential expression for Trem2-high vs Trem2-low
DEG_Trem2 <- FindMarkers(
  Macs,
  ident.1 = "Trem2-high",
  ident.2 = "Trem2-low",
  verbose = TRUE,
  logfc.threshold = 0.05
)

# Add gene names as a column
DEG_Trem2$gene <- rownames(DEG_Trem2)

# Export to CSV
write.csv(DEG_Trem2, "DEG_Trem2_wSymbols.csv", row.names = FALSE)

# Add -log10 adjusted p-value for plotting
DEG_Trem2$neg_adj.pVal <- -log10(DEG_Trem2$p_val_adj)

# Define threshold for significance
pThresh <- 2  # equivalent to adjusted p < 0.01

# Filter out Trem2 from the volcano plot (but still in CSV)
DEG_Trem2_plot <- DEG_Trem2[DEG_Trem2$gene != "Trem2", ]

# Base R volcano plot
plot(
  DEG_Trem2_plot$avg_log2FC,
  DEG_Trem2_plot$neg_adj.pVal,
  col = ifelse(DEG_Trem2_plot$neg_adj.pVal > pThresh, "red3", "black"),
  ylim = c(0, 50),
  pch = ifelse(DEG_Trem2_plot$neg_adj.pVal > pThresh, 19, 1),
  cex = 1.2,
  xlab = "avg_log2FC",
  ylab = "-log10(adj.p-value)",
  main = "Volcano Plot (Trem2 DEGs)"
)

abline(v = seq(-6, 6, by = 2), col = "lightgray")
abline(h = seq(0, 50, by = 10), col = "lightgray")
abline(h = pThresh, col = "black", lty = 2)

# Highlight up/down DEGs
points(
  DEG_Trem2_plot$avg_log2FC,
  DEG_Trem2_plot$neg_adj.pVal,
  col = ifelse(DEG_Trem2_plot$neg_adj.pVal > pThresh,
               ifelse(DEG_Trem2_plot$avg_log2FC > 0, "red3", "dodgerblue3"),
               "black"),
  pch = ifelse(DEG_Trem2_plot$neg_adj.pVal > pThresh, 19, 1),
  cex = ifelse(DEG_Trem2_plot$neg_adj.pVal > pThresh, 1.4, 1.2)
)

# ggplotly Volcano Plot (interactive)
library(ggplot2)
library(plotly)

ggplotly(
  ggplot(
    data = DEG_Trem2_plot,
    aes(
      x = avg_log2FC,
      y = neg_adj.pVal,
      text = gene,
      colour = ifelse(neg_adj.pVal > pThresh, "darkblue", "darkred")
    )
  ) +
    geom_point(shape = 21, size = 3, fill = "white") +
    geom_hline(yintercept = 0) +
    theme_minimal() +
    labs(title = "Volcano Plot (without Trem2)", x = "avg_log2FC", y = "-log10 adj.p")
)

### LPM vs SPM
LPM_markers <- c("Adgre1", "Timd4", "Gata6", "Icam2")
SPM_markers <- c("Ly6c2", "Ccr2", "Cd14", "Fcgr1")

Myeloids <- AddModuleScore(
  Myeloids,
  features = list(LPM_markers),
  name = "LPM_score"
)

Myeloids <- AddModuleScore(
  Myeloids,
  features = list(SPM_markers),
  name = "SPM_score"
)

FeaturePlot(Myeloids, features = c("LPM_score1", "SPM_score1"), reduction = "umap")

Myeloids$Macrophage_Type <- ifelse(
  Myeloids$LPM_score1 > Myeloids$SPM_score1,
  "LPM",
  "SPM"
)

library(dplyr)

prop_df <- Myeloids@meta.data %>%
  group_by(Combined.HTO_group, Macrophage_Type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

ggplot(prop_df, aes(x = Combined.HTO_group, y = freq, fill = Macrophage_Type)) +
  geom_bar(stat = "identity")

VlnPlot(
  Myeloids,
  features = c("Tnf", "Il1b", "Il10", "Tgfb1"),
  group.by = "Macrophage_Type"
)
Idents(Myeloids) <- Myeloids$Macrophage_Type
ligands <- c("Tnf", "Il1b", "Il10", "Tgfb1")

Myeloids <- PrepSCTFindMarkers(Myeloids)

DE_LPM_vs_SPM <- FindMarkers(
  Myeloids,
  ident.1 = "SPM",
  ident.2 = "LPM",
  features = ligands,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0
)

DE_LPM_vs_SPM$gene <- rownames(DE_LPM_vs_SPM)

DE_LPM_vs_SPM <- DE_LPM_vs_SPM %>%
  dplyr::select(
    gene,
    avg_log2FC,
    p_val,
    p_val_adj,
    pct.1,
    pct.2
  )

write.csv(
  DE_LPM_vs_SPM,
  file = "SPM_vs_LPM_cytokine_expression.csv",
  row.names = FALSE
)

## DEG analysis for LPM: Ascites vs OmMet
# 1. Subset the object to only include LPMs
LPM_subset <- subset(Myeloids, subset = Macrophage_Type == "LPM")

# 2. Set identities to the HTO group
Idents(LPM_subset) <- "Combined.HTO_group"

# 3. Run FindMarkers comparing the two specific HTOs
# Adjust ident.1 and ident.2 based on which one you want as the 'test' vs 'control'
LPM_subset <- SCTransform(LPM_subset, assay = "RNA", verbose = FALSE)

# 2. Now run PrepSCT and FindMarkers
LPM_subset <- PrepSCTFindMarkers(LPM_subset)

DEG_HTO_LPM <- FindMarkers(
  LPM_subset, 
  ident.1 = "asc:1_AscMet", 
  ident.2 = "om:2_OmMet", 
  assay = "SCT", # Explicitly use the unified SCT assay
  recorrect_erythrocytes = FALSE # Turn off if not needed for the subset
)

# 4. Add the negative log10 p-value for plotting
DEG_HTO_LPM$neg_log10_padj <- -log10(DEG_HTO_LPM$p_val_adj)
DEG_HTO_LPM$gene <- rownames(DEG_HTO_LPM)

# Categorize for coloring
DEG_HTO_LPM <- DEG_HTO_LPM %>%
  mutate(diffexpressed = case_when(
    avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up",
    avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down",
    TRUE ~ "No"
  ))

# Select top genes for labeling (15 up, 15 down)
top_labels <- rbind(
  DEG_HTO_LPM %>% filter(diffexpressed == "Up") %>% slice_min(p_val_adj, n = 25),
  DEG_HTO_LPM %>% filter(diffexpressed == "Down") %>% slice_min(p_val_adj, n = 15)
)

# Create Plot
ggplot(DEG_HTO_LPM, aes(x = avg_log2FC, y = neg_log10_padj, fill = diffexpressed)) +
  geom_point(aes(size = neg_log10_padj), shape = 21, color = "black", stroke = 0.3, alpha = 0.8) +
  
  # Aesthetic Theme
  scale_fill_manual(values = c("Down" = "black", "No" = "grey80", "Up" = "firebrick3")) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  theme_bw(base_size = 14) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", size = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.4) +
  
  # Labeling
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5,
    max.overlaps = 30,
    segment.size = 0.2
  ) +
  
  # Labels
  labs(
    title = "LPMs: AscMet vs. OmMet",
    subtitle = "Positive FC = Enriched in AscMet",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P-value),
    fill = "Significance"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
    axis.text = element_text(color = "black")
  )

write.csv(DEG_HTO_LPM, file = "DEG_AscMet_vs_OmMet_LPM.csv", row.names = FALSE)

### EXPORT SOURCE DATA
library(tidyverse)
library(writexl)

# ==========================================================
# Figure 5D: Broad Cell ID Proportions
# ==========================================================
# Convert your prop.table into a clean dataframe
fig5d_broad_props <- as.data.frame(
  prop.table(table(Ova_merged_SCT@meta.data$Broad_cell.ID, Ova_merged_SCT@meta.data$Combined.HTO_group.order), margin = 2)
) %>%
  rename(Broad_Cell_Type = Var1, Experimental_Group = Var2, Proportion = Freq) %>%
  # Filter out 0s and format neatly
  filter(Proportion > 0) %>%
  arrange(Experimental_Group, desc(Proportion))

# ==========================================================
# Figure 5E & 5F: Myeloid UMAPs and LPM/SPM Feature Scores
# ==========================================================
# Extract the UMAP coordinates and the exact module scores used for the feature plots
fig5e_5f_myeloid_umap <- FetchData(
  Myeloids, 
  vars = c("UMAP_1", "UMAP_2", "Combined.HTO_group.order", "LPM_score1", "SPM_score1", "Macrophage_Type")
) %>%
  rownames_to_column(var = "Cell_Barcode") %>%
  rename(Experimental_Group = Combined.HTO_group.order) %>%
  arrange(Experimental_Group)

# ==========================================================
# Figure 5G: LPM / SPM Ratios
# ==========================================================
# Grab the prop_df you already calculated
fig5g_lpm_spm_ratios <- prop_df %>%
  rename(
    Experimental_Group = Combined.HTO_group,
    Cell_Count = n,
    Proportion = freq
  ) %>%
  arrange(Experimental_Group, Macrophage_Type)

# ==========================================================
# Figure 5I: LPM AscMet vs OmMet Volcano Plot
# ==========================================================
# Clean up the DEG table to provide the exact coordinates and p-values
fig5i_lpm_volcano <- DEG_HTO_LPM %>%
  select(gene, p_val, avg_log2FC, p_val_adj, neg_log10_padj, diffexpressed) %>%
  arrange(p_val_adj)

# ==========================================================
# Figure 5L: T/NK Cell Proportions
# ==========================================================
# Convert the T.NK prop.table into a clean dataframe
fig5l_tnk_props <- as.data.frame(
  prop.table(table(T.NK@meta.data$seurat_clusters, T.NK@meta.data$Combined.HTO_group.order), margin = 2)
) %>%
  rename(T_NK_Cluster = Var1, Experimental_Group = Var2, Proportion = Freq) %>%
  filter(Proportion > 0) %>%
  arrange(Experimental_Group, T_NK_Cluster)

# ==========================================================
# Package and Export
# ==========================================================
figure5_source_data <- list(
  "Fig5D_Broad_Cell_Props" = fig5d_broad_props,
  "Fig5E_5F_Myeloid_UMAP"  = fig5e_5f_myeloid_umap,
  "Fig5G_LPM_SPM_Ratios"   = fig5g_lpm_spm_ratios,
  "Fig5I_LPM_Volcano"      = fig5i_lpm_volcano,
  "Fig5L_TNK_Proportions"  = fig5l_tnk_props
)

write_xlsx(figure5_source_data, "Source_Data_Figure5_Immune_Microenvironment.xlsx")

print("Figure 5 Source Data exported successfully!")
