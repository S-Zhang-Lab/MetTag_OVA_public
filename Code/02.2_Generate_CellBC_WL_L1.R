# ==== Load Cellranger output, perform QC and write Cell BC =====
# Step 1 -load 10x output. "filtered_feature_bc_matrix" is under 
# "outs" folder from Cellranger output folder. 
# bioHCP path is: /endosome/archive/pathology/SiZhang_lab/shared/UTSW_CITE_RAW/R02_EA/cellranger.outs/R02_EA_L1_LARRYreseq/outs
# This code is for L1. Need to change path for L2
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

base_dir <- "C:/Users/s208205/Downloads/R_projects/MetTag_BC_OVA"
figures_dir <- file.path(base_dir, "L1_Figures")

MetTag.ova <- Read10X(
  data.dir = file.path(base_dir, "/DATA/R02_EA_L1_LARRYreseq/filtered_feature_bc_matrix"))
str(MetTag.ova) # there should be two assays: Gene Expression, Antibody, and Custom (LARRY BC)

# establish Seurat project
L1 <- CreateSeuratObject(counts = MetTag.ova$"Gene Expression", project = "L1_Om")
L1$CITE <- CreateAssayObject(counts = MetTag.ova$"Antibody Capture")

# ADT/HTO processing. Remove ADT columne if there was ADT antibody used. 
L1.CITE <- t(as.data.frame(L1$"CITE"@counts))
colnames(L1.CITE)

L1.HTOs <- L1.CITE[, c(47:49)]
L1[["HTO"]] <- CreateAssayObject(counts = t(L1.HTOs))

# Perform QC
L1 <- PercentageFeatureSet(L1,pattern = "^mt-", col.name = "percent.mt")
meta <- L1@meta.data

# Visualize QC metrics as a violin plot
pdf(file.path(figures_dir, "QC1.pdf"),width=6,height=10)
VlnPlot(L1, features = c(
  "nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_HTO"), ncol = 2)
dev.off()

pdf(file.path(figures_dir, "QC2.pdf"),width=4,height=4)
VlnPlot(L1, features = c(
  "nFeature_HTO", "nCount_HTO"), ncol = 2) # 
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf(file.path(figures_dir, "QC_Featureplot_RNA.pdf"),width=10,height=4,paper='special')
plot1 <- FeatureScatter(L1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(L1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

pdf(file.path(figures_dir, "QC_Featureplot_HTO_BC.pdf"),width=6,height=6,paper='special')
plot3 <- FeatureScatter(L1, feature1 = "nFeature_HTO", feature2 = "nCount_HTO")
plot3
dev.off()

# Add ratio calculation between nCount and nFeature. 
L1@meta.data$RNA_ratio <- L1@meta.data$nCount_RNA/L1@meta.data$nFeature_RNA
L1@meta.data$HTO_ratio <- L1@meta.data$nCount_HTO/L1@meta.data$nFeature_HTO

# Removing unwanted cells from the data based on above QC plots
L1 <- subset(L1, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 50)

# Normalize HTOs and subset singlets
L1 <- NormalizeData(L1, assay = "HTO", normalization.method = "CLR")
L1 <- HTODemux(L1, assay = "HTO", positive.quantile = 0.99)
meta <- L1@meta.data
table(meta$HTO_classification.global)

# Group cells based on the max HTO signal
DefaultAssay(L1) <- "HTO"
Idents(L1) <- "HTO_maxID"
pdf(file.path(figures_dir, "QC_RidgePlot_HTO3-5.pdf"),width=10,height=8,paper='special')
RidgePlot(L1, assay = "HTO", features = rownames(L1[["HTO"]])[1:4], ncol = 2)
dev.off()

pdf(file.path(figures_dir, "QC_Scatter_HTO.pdf"),width=10,height=8,paper='special')
plot1 <- FeatureScatter(L1, feature1 = "M-HTO-1", feature2 = "M-HTO-2")
plot2 <- FeatureScatter(L1, feature1 = "M-HTO-3", feature2 = "M-HTO-4")
plot3 <- FeatureScatter(L1, feature1 = "M-HTO-1", feature2 = "M-HTO-3")
plot4 <- FeatureScatter(L1, feature1 = "M-HTO-2", feature2 = "M-HTO-4")
plot1 + plot2 + plot3 + plot4
dev.off()

# subset singlets
Idents(L1) <- "HTO_classification.global"
L1.singlet <- subset(L1, idents = "Singlet")
meta.L1.singlet <- L1.singlet@meta.data

L1.CellBC_WL <- rownames(L1.singlet@meta.data)
L1.CellBC_WL <- sub("-1$", "", L1.CellBC_WL) # Remove the trailing "-1"

writeLines(L1.CellBC_WL, con = "./DATA/BioHPC.files/L1_CellBC_WL.txt")  # this is the file for 03_correctBC_parallel_v3.py

saveRDS(L1.singlet, file = "./DATA/RDS/L1.singlet.rds")

