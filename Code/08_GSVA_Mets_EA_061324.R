############################################################## 
#          GSVA analysis for R02 CITE-seq                    # 
############################################################## 
# loading packages and gene sets 
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(GSEABase) 
library(GSVA) 
library(GSVAdata) 
library(Biobase) 
library(limma) 
library(RColorBrewer) 
library(parallel) 

setwd("C:/Users/Zhang Lab/Desktop/Emma/R01_and_2_w_LARRY/GSVA")
Mets <- readRDS("Mets_SCT_EA_062024.rds")

# Find tumor cluster markers 
DefaultAssay(Mets)
Idents(Mets) <- "seurat_clusters"
Mets <- PrepSCTFindMarkers(Mets)

Mets.markers <- FindAllMarkers(Mets, only.pos = TRUE)
Mets.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

Mets.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 50) %>%
  ungroup() -> top50

write.csv(top50, "Met_cluster_markerTop50.csv")

# get gene sets
gs=getGmt("msigdb.v2022.1.Mm.symbols.gmt") 

# v7 GSEA MSigDB gene sets, gene symbols, All genesets 
summary(gs) 
gs.full <- gs 

# run ssGSEA on all genesets, if not previously run. Otherwise can load directly. See the next step.  
ssGSEA_result_gs.full = gsva(as.matrix(Mets@assays$SCT$data), gs.full, min.sz=5, max.sz=500, verbose=TRUE, method="ssgsea") 
write.csv(ssGSEA_result_gs.full, file = "R02_ssGSEA_result_gsFull_integrated_Mets.csv") 
ssGSEA_result_gs.full <- read.csv("./08_GSVA_Met_Figures/R02_ssGSEA_result_gsFull_integrated_Mets.csv", row.names = 1)

# Min-Max scaling 
# Load necessary library
library(scales)

# Assuming ssGSEA_result_gs.full contains GSVA scores
# Apply Min-Max scaling to transform scores to range [0, 1]
min_max_scaled_scores <- apply(ssGSEA_result_gs.full, 2, function(x) {
  rescaled_x <- scales::rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE))
  return(rescaled_x)
})


# Create a new Seurat object with the Min-Max scaled GSVA scores
gsvasc <- CreateSeuratObject(counts = min_max_scaled_scores, project = "ssGSEA_result")
gsvasc 

any(gsvasc@assays$RNA$counts < 0) # check whether data contain negative values
any(gsvasc@assays$RNA$counts == 0)

gsva.meta <- gsvasc@meta.data

# Merge the meta
gsvasc@meta.data <- Mets@meta.data

head(gsvasc@meta.data) 
colnames(gsvasc@meta.data)

# Normalization
gsvasc <- NormalizeData(gsvasc)

# find variable features
gsvasc <- FindVariableFeatures(gsvasc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(gsvasc)
gsvasc <- ScaleData(gsvasc, features = all.genes)

plot1 <- VariableFeaturePlot(gsvasc)
plot1

# Set seed and create copy of R object
set.seed(42)
gsvasc.2 <- gsvasc

# Examine top marker pathways for each met cluster
Idents(gsvasc) <- "seurat_clusters"
gs.cluster.markers <- FindAllMarkers(gsvasc, only.pos = TRUE) 
gs.cluster.markers  %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> gs.cluster.markers.top50 
write.csv(as.matrix(gs.cluster.markers.top50), file = "gs.cluster_markers_top50.csv") 

# check individual pathway expression
Idents(gsvasc) <- "Combined.HTO_group"
VlnPlot(gsvasc, features = "RAFFEL-VEGFA-TARGETS-DN")
VlnPlot(gsvasc, features = "GOCC-HAPTOGLOBIN-HEMOGLOBIN-COMPLEX")
VlnPlot(gsvasc, features = "GOBP-MONOCYTE-CHEMOTACTIC-PROTEIN-1-PRODUCTION")
VlnPlot(gsvasc, features = "GOMF-EXTRACELLULAR-MATRIX-STRUCTURAL-CONSTITUENT")
VlnPlot(gsvasc, features = "GOBP-ARACHIDONIC-ACID-SECRETION", split.by = "Combined.HTO_group")
VlnPlot(gsvasc, features = c("REACTOME-INITIAL-TRIGGERING-OF-COMPLEMENT", "GOBP-NEGATIVE-REGULATION-OF-CD8-POSITIVE-ALPHA-BETA-T-CELL-ACTIVATION", "GOBP-COMPLEMENT-RECEPTOR-MEDIATED-SIGNALING-PATHWAY"))

# save progress
saveRDS(gsvasc, file = "08_Mets_gsvasc.rds")

# ========== Global DE pathway Identify the DE pathways between OmMet and AscMet  =============

Idents(gsvasc) <- "Combined.HTO_group" 
table(Idents(gsvasc))
DEG_OmMet_AscMet_gsvasc <- FindMarkers(gsvasc, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", verbose = TRUE, logfc.threshold = 0.05)

# ========== Global DE pathway Identify the DE pathways between Expanded and Non-expanded  =============

Idents(gsvasc) <- "Combined.HTO_group" 
gsvasc.Om <- subset(gsvasc, idents = "om:2_OmMet")
Om.meta <- gsvasc.Om@meta.data
barcode_counts <- table(Om.meta$Combined.LARRY_group)
Om.meta$Expansion_value <- "Mid"
expanded_barcodes <- names(barcode_counts[barcode_counts >= 10])
non_expanded_barcodes <- names(barcode_counts[barcode_counts == 1])
not_detected_barcodes <- "LARRY_ratio_less 2-om"
Om.meta$Expansion_value[Om.meta$Combined.LARRY_group %in% expanded_barcodes] <- "Expanded"
Om.meta$Expansion_value[Om.meta$Combined.LARRY_group %in% non_expanded_barcodes] <- "Non-expanded"
Om.meta$Expansion_value[Om.meta$Combined.LARRY_group %in% not_detected_barcodes] <- "Not-detected"

# Add the updated metadata back to the Seurat object
gsvasc.Om <- AddMetaData(object = gsvasc.Om, metadata = Om.meta)
# Verify the new metadata slot
table(gsvasc.Om@meta.data$Expansion_value)
Idents(gsvasc.Om) <- "Expansion_value"
DEG_OmMet_expansion <- FindMarkers(gsvasc.Om, ident.1 = "Expanded", ident.2 = "Non-expanded", verbose = TRUE, logfc.threshold = 0.05)

# save both comparisons
write.csv(DEG_OmMet_expansion, file = "DEG_OmMet_expansion_pathways.csv")
write.csv(DEG_OmMet_AscMet_gsvasc, file = "08_DEG_OmMet_AscMet_gsvasc.csv")

# OmMet vs AscMet contains a lot of pathways, filter out most significant ones for plotting
DEG_OmMet_AscMet_gsvasc_Top200_sorted <- DEG_OmMet_AscMet_gsvasc %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 200) 
write.csv(DEG_OmMet_AscMet_gsvasc_Top200_sorted, file = "DEG_OmMet_AscMet_gsvasc_Top200_sorted.csv") 

# save workspace
save.image("08_GSVA_Mets.RData")

###Plot GSVA data in barplot/cleveland plot style
#1. combine DEGS lists acquired for control & ABX between different time points
#2. then filter out immune-related gene sets
#3. then create barplot or cleveland plot show immune-related GS enrichment for each condition between time points

library(ggplot2)
library(dplyr)
library(readxl)

#combine DEGS lists acquired for control & ABX between different time points + expanded vs. non expanded
#import excel files generated when measuring significant DEGS between time points per conditions 
#(removed all gene sets with adj p value > 0.05)
DEG_AscMet <- read_excel("DEG_AscMet_gsvasc_sorted.xlsx")
DEG_OmMet <- read_excel("DEG_OmMet_gsvasc_sorted.xlsx")

DEG_count1 <- read_excel("DEG_OmMet_expansion_count1.xlsx")
DEG_Expanded <- read_excel("DEG_OmMet_expansion_top3.xlsx")

#change column name from "...1" to "gene_set"
#the "...1" column name only happens if you import a DEG table that you had exported after running FindMarkers
colnames(DEG_count1)[colnames(DEG_count1)=="...1"]<-"gene_set"
colnames(DEG_Expanded)[colnames(DEG_Expanded)=="...1"]<-"gene_set"

colnames(DEG_OmMet)[colnames(DEG_OmMet)=="...1"]<-"gene_set"
colnames(DEG_AscMet)[colnames(DEG_AscMet)=="...1"]<-"gene_set"

#add condition information to each table
DEG_count1$Condition <- "01_Non-expanded"
DEG_Expanded$Condition <- "02_Expanded"

DEG_AscMet$Condition <- "01_AscMet"
DEG_OmMet$Condition <- "02_OmMet"

#combine tables
comb_DEGS <- rbind(DEG_count1, DEG_Expanded)

comb_DEGS.tx <- rbind(DEG_AscMet, DEG_OmMet)

#create Cleveland/Lollipop plot with combined DEGS information
#reduce list for plot by filtering out DEGS with > abs value 0.25 avg_log2FC
comb_DEGS.tx <- comb_DEGS.tx %>% filter(abs(avg_log2FC) >= 0.25)

#reduce list for plot by filtering out DEGS with > abs value 0.1 avg_log2FC
# for expansion comparison
comb_DEGS <- comb_DEGS %>% filter(abs(avg_log2FC) >= 0.1)

# Define colors for each condition
condition_colors <- c("01_Non-expanded" = "gray", "02_Expanded" = "purple")
condition_colors.tx <- c("01_AscMet" = "blue", "02_OmMet" = "red")

# Create the plot with specific colors for each condition
#you can edit this part of the code if you'd rather plot p val on the x-axis rather than Log FC
#alternatively, you could have Log FC on the x-axis but have the bars colored by p value (I tried this but thought it looked weird/not very infomative because all the p-values were very low for my filtered gene set list)
ggplot(comb_DEGS, aes(x = avg_log2FC, y = reorder(gene_set, -avg_log2FC))) +
  geom_segment(aes(xend = 0, yend = reorder(gene_set, avg_log2FC), color = Condition), linewidth = 0.5, size = 2.5) +
  geom_point(aes(color = Condition)) +
  scale_color_manual(values = condition_colors) +  # Assign specific colors for each condition
  facet_grid(Condition ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Average Log2FC", y = "Gene Set", color = "Condition") +
  theme_minimal()
# export 4x16 pdf

ggplot(comb_DEGS.tx, aes(x = avg_log2FC, y = reorder(gene_set, -avg_log2FC))) +
  geom_segment(aes(xend = 0, yend = reorder(gene_set, avg_log2FC), color = Condition), linewidth = 0.5, size = 2.5) +
  geom_point(aes(color = Condition)) +
  scale_color_manual(values = condition_colors.tx) +  # Assign specific colors for each condition
  facet_grid(Condition ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Average Log2FC", y = "Gene Set", color = "Condition") +
  theme_minimal()
# export 8x16 pdf

### GSEA on Gbp2b high vs low cells

# get gene sets - GO terms 
gs=getGmt("m5.go.v2024.1.Mm.symbols.gmt") 

# v7 GSEA MSigDB gene sets, gene symbols, All genesets 
summary(gs) 
gs.full <- gs 
# run ssGSEA on all genesets, if not previously run. Otherwise can load directly. See the next step.  
ssGSEA_result_gs.full = gsva(as.matrix(Mets@assays$SCT$data), gs.full, min.sz=5, max.sz=500, verbose=TRUE, method="ssgsea") 
write.csv(ssGSEA_result_gs.full, file = "R02_ssGSEA_result_gsGO_integrated_Mets.csv") 

Idents(gsvasc) <- "Gbp2b_status"
DEG_gsvasc_Gbp2bhigh <- FindMarkers(gsvasc, ident.1 = "Gbp2b(high)", ident.2 = "Gbp2b(low)", verbose = TRUE, logfc.threshold = 0.05)
write.csv(as.matrix(DEG_gsvasc_Gbp2bhigh), file = "DEG_gsvasc_Gbp2bHigh_vs_Low.csv")

