### Step7.3.1 analysis UTSW02 - cell annotation ###

# Load Packages + DATA 
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

#Load in data and examine structure (note:Broad_cell.ID added in step 07.3.2)
# for OG object before these examining cell types:
# Ova_merged_SCT <- readRDS("Step7.1_Ova_merged_SCT_annotated.rds")

Ova_merged_SCT <- readRDS("Ova_merged_SCT_EA_061824.rds")
Idents(Ova_merged_SCT) <- "seurat_clusters"
colnames(Ova_merged_SCT@meta.data)
meta_SCT <- Ova_merged_SCT@meta.data
meta_SCT$seurat_clusters %>% table()
meta_SCT$Combined.HTO_group %>% table()
meta_SCT$Broad_cell.ID %>% table()
DimPlot(Ova_merged_SCT)
DimPlot(Ova_merged_SCT, split.by = "Combined.HTO_group")

### 1_CELL IDENTIFICATION/GENERAL ###
#Examine Top Expressed Marker Genes Among SCT Clusters
##Note: previously used resolution yields 42 clusters (a lot) but we mostly focused on Broad_cell.ID analysis
Idents(Ova_merged_SCT) <- "seurat_clusters"
Ova_merged_SCT <- PrepSCTFindMarkers(Ova_merged_SCT)
pbmc.markers <- FindAllMarkers(Ova_merged_SCT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, file = "Ova_merged_SCT_RNA_markergenes.csv")
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Ova_merged_SCT, features = top10$gene) + NoLegend()+scale_fill_gradientn(colors = c("blue", "white", "red")) 
# save pdf 26x24
# visualize immune cells (Ptprc+ Krt8- Krt18-), cancer cells and other stroma (Ptprc- Krt8+, Krt18+), cancer cells only (Ptprc- Egfp+ Krt8+ Krt18+)
FeaturePlot(Ova_merged_SCT, features = c("Ptprc", "Krt8","Krt18", "EGFP"), min.cutoff = "q05", max.cutoff = "q95")
#save pdf 6x7

# Note: the following violin plot were used in conjunction with the top10 gene heatmap/marker genes csv + PanglaoDB
# to assign cell type identities to clusters

# ID efferocytotic macs - clusters 5 and 8
VlnPlot(Ova_merged_SCT, features = c("Trem2", "Ccr5", "C1qa", "Il18bp", "Adgre1", "Marco"))

# ID T cells 
VlnPlot(Ova_merged_SCT, features = c("Cd3e", "Cd8a", "Cd4", "Foxp3", "Cd44", "Pdcd1"))

# ID NK and NKT cells - clusters 21,23
VlnPlot(Ova_merged_SCT, features = c("Itga2", "Klrb1c", "Ncam1", "Gzma", "Cd44", "Cd3e"))

#ID gdT cells - clusters 30,31
VlnPlot(Ova_merged_SCT, features = c("Tcrg-C1", "Trgv2", "Tcrg-C2", "Blk", "Rora"))
# save pdf 7x14

# B cells cluster 2, 3 and more
VlnPlot(Ova_merged_SCT, features = c("Fcer2a", "Ighd", "Bhlhe41", "Zbtb32", "Ighm", "Vpreb3"))

# Pericytes
VlnPlot(Ova_merged_SCT, features = c("Pdgfrb", "Acta2", "Cspg4", "Mcam", "Angpt1", "Des"))

# Pdgfrb fibroblasts
VlnPlot(Ova_merged_SCT, features = c("Acta2", "Kcnk3", "Pdgfrb", "Higd1b", "Cox4i2", "Notch3", "Postn"))
# save pdf 10x14
# macs vs. mono
VlnPlot(Ova_merged_SCT, features = c("Csf1r", "Itgam", "Ccr5", "Adgre1", "Itgax", "Ccr2", "Ly6c1"))

#Make copies of original file 
Ova_merged_SCT.2 <- Ova_merged_SCT

#Compare Cluster Proportions Between Control and Experimental 
table(Ova_merged_SCT@meta.data$seurat_clusters,Ova_merged_SCT@meta.data$Combined.HTO_group)
freq_table <- prop.table(x = table(Ova_merged_SCT@meta.data$seurat_clusters,Ova_merged_SCT@meta.data$Combined.HTO_group),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(Ova_merged_SCT@meta.data$seurat_clusters) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)
# export one at 7x10 to see x axis labels correctly
# export one 7x5 - more aesthetically pleasing
# Ifng response accross groups

Idents(Ova_merged_SCT) <- "Broad_cell.ID"
VlnPlot(Ova_merged_SCT, features = c("Stat1", "Irf1", "Ifngr1", "Gbp2b", "Ifng", "Marco", "Ms4a8a", "Cd300lb", "Slfn1"))

# NOTE 2: SZ conducted cancer cell subsetting > do not use EA subsetting in 7.3.2
