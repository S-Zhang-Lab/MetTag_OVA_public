############################################################## 
#          GSVA analysis for R02 CITE-seq - Part2            # 
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

#set working directory
setwd("C:/Users/Zhang Lab/Desktop/Emma/R01_and_2_w_LARRY/Step8_GSVA_DATA")
Mets <- readRDS("Mets_SCT_EA_062024.rds")

# get gene sets
# Run individual Hallmark, Reactome, GO sets instead of everything together
gs.h=getGmt("mh.all.v2024.1.Mm.symbols (1).gmt") 
gs.r=getGmt("m2.cp.reactome.v2024.1.Mm.symbols.gmt")
gs.go=getGmt("m5.go.v2024.1.Mm.symbols.gmt")
summary(gs.h) 
summary(gs.r)
summary(gs.go)

# run ssGSEA on Hallmark gene set  
ssGSEA_result_gs.h = gsva(as.matrix(Mets@assays$SCT$data), gs.h, min.sz=5, max.sz=500, verbose=TRUE, method="ssgsea") 
write.csv(ssGSEA_result_gs.h, file = "R02_ssGSEA_result_gsH_integrated_Mets.csv") 
ssGSEA_result_gs.h <- read.csv("R02_ssGSEA_result_gsH_integrated_Mets.csv", row.names = 1)

# run ssGSEA on GO gene set  
ssGSEA_result_gs.go = gsva(as.matrix(Mets@assays$SCT$data), gs.go, min.sz=5, max.sz=500, verbose=TRUE, method="ssgsea") 
write.csv(ssGSEA_result_gs.go, file = "R02_ssGSEA_result_gsGO_integrated_Mets.csv") 

ssGSEA_result_gs.full <- ssGSEA_result_gs.go

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
write.csv(as.matrix(gs.cluster.markers.top50), file = "gs.cluster_GO_markers_top50.csv") 

# check individual pathway expression
Idents(gsvasc) <- "Combined.HTO_group"
VlnPlot(gsvasc, features = "RAFFEL-VEGFA-TARGETS-DN")
VlnPlot(gsvasc, features = "GOCC-HAPTOGLOBIN-HEMOGLOBIN-COMPLEX")
VlnPlot(gsvasc, features = "GOBP-MONOCYTE-CHEMOTACTIC-PROTEIN-1-PRODUCTION")
VlnPlot(gsvasc, features = "GOMF-EXTRACELLULAR-MATRIX-STRUCTURAL-CONSTITUENT")
VlnPlot(gsvasc, features = "GOBP-ARACHIDONIC-ACID-SECRETION", split.by = "Combined.HTO_group")
VlnPlot(gsvasc, features = c("REACTOME-INITIAL-TRIGGERING-OF-COMPLEMENT", "GOBP-NEGATIVE-REGULATION-OF-CD8-POSITIVE-ALPHA-BETA-T-CELL-ACTIVATION", "GOBP-COMPLEMENT-RECEPTOR-MEDIATED-SIGNALING-PATHWAY"))

# save progress
saveRDS(gsvasc, file = "09_Mets_gsvasc_GO.rds")

# ========== Global DE pathway Identify the DE pathways between OmMet and AscMet  =============

Idents(gsvasc) <- "Combined.HTO_group" 
table(Idents(gsvasc))
DEG_OmMet_AscMet_gsvasc <- FindMarkers(gsvasc, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_OmMet_AscMet_gsvasc, file = "09_GO_OmMet_AscMet_gsvasc.csv")

###Plot GSVA data in barplot/cleveland plot style
#1. combine DEGS lists acquired for control & ABX between different time points
#2. then filter out immune-related gene sets
#3. then create barplot or cleveland plot show immune-related GS enrichment for each condition between time points

library(ggplot2)
library(dplyr)
library(readxl)

#combine DEGS lists acquired for AscMet and OmMet
#import excel files generated when measuring significant DEGS between time points per conditions 
#(removed all gene sets with adj p value > 0.05, log2FC <0.25)
DEG_AscMet <- read_excel("GOBP_Up_AscMet.xlsx")
DEG_OmMet <- read_excel("GOBP_up_OmMet.xlsx")
colnames(DEG_OmMet)[colnames(DEG_OmMet)=="...1"]<-"gene_set"
colnames(DEG_AscMet)[colnames(DEG_AscMet)=="...1"]<-"gene_set"
DEG_AscMet$Condition <- "01_AscMet"
DEG_OmMet$Condition <- "02_OmMet"
comb_DEGS.tx <- rbind(DEG_AscMet, DEG_OmMet)

#create Cleveland/Lollipop plot with combined DEGS information
#reduce list for plot by filtering out DEGS with > abs value 0.25 avg_log2FC
comb_DEGS.tx <- comb_DEGS.tx %>% filter(abs(avg_log2FC) >= 0.3)
condition_colors.tx <- c("01_AscMet" = "blue", "02_OmMet" = "red")
# Create the plot with specific colors for each condition
#you can edit this part of the code if you'd rather plot p val on the x-axis rather than Log FC
#alternatively, you could have Log FC on the x-axis but have the bars colored by p value (I tried this but thought it looked weird/not very infomative because all the p-values were very low for my filtered gene set list)
ggplot(comb_DEGS.tx, aes(x = avg_log2FC, y = reorder(gene_set, -avg_log2FC))) +
  geom_segment(aes(xend = 0, yend = reorder(gene_set, avg_log2FC), color = Condition), linewidth = 0.5, size = 2.5) +
  geom_point(aes(color = Condition)) +
  scale_color_manual(values = condition_colors.tx) +  # Assign specific colors for each condition
  facet_grid(Condition ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Average Log2FC", y = "Gene Set", color = "Condition") +
  theme_minimal()
# export 8x16 pdf


