# ==== Load Cellranger output, perform QC and write Cell BC =====
# Step 1 -load 10x output. "filtered_feature_bc_matrix" is under 
# "outs" folder from Cellranger output folder. 
# bioHCP path is: /endosome/archive/pathology/SiZhang_lab/shared/UTSW_CITE_RAW/R02_EA/cellranger.outs/R02_EA_L2_LARRYreseq/outs
# This code is for L2. Need to change path for L2
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

base_dir <- "C:/Users/s208205/Downloads/R_projects/MetTag_BC_OVA"
figures_dir <- file.path(base_dir, "L2_Figures")

MetTag.ova <- Read10X(
  data.dir = file.path(base_dir, "/DATA/R02_EA_L2_LARRYreseq/filtered_feature_bc_matrix"))
str(MetTag.ova) # there should be two assays: Gene Expression, Antibody, and Custom (LARRY BC)

# establish Seurat project
L2 <- CreateSeuratObject(counts = MetTag.ova$"Gene Expression", project = "L2_Asc")
L2$CITE <- CreateAssayObject(counts = MetTag.ova$"Antibody Capture")

# ADT/HTO processing. Remove ADT columne if there was ADT antibody used. 
L2.CITE <- t(as.data.frame(L2$"CITE"@counts))
colnames(L2.CITE)

L2.HTOs <- L2.CITE[, c(45:52)]
L2[["HTO"]] <- CreateAssayObject(counts = t(L2.HTOs))

# Perform QC
L2 <- PercentageFeatureSet(L2,pattern = "^mt-", col.name = "percent.mt")
meta <- L2@meta.data

# Visualize QC metrics as a violin plot
pdf(file.path(figures_dir, "QC1.pdf"),width=6,height=10)
VlnPlot(L2, features = c(
  "nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_HTO"), ncol = 2)
dev.off()

pdf(file.path(figures_dir, "QC2.pdf"),width=4,height=4)
VlnPlot(L2, features = c(
  "nFeature_HTO", "nCount_HTO"), ncol = 2) # 
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf(file.path(figures_dir, "QC_Featureplot_RNA.pdf"),width=10,height=4,paper='special')
plot1 <- FeatureScatter(L2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(L2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

pdf(file.path(figures_dir, "QC_Featureplot_HTO_BC.pdf"),width=6,height=6,paper='special')
plot3 <- FeatureScatter(L2, feature1 = "nFeature_HTO", feature2 = "nCount_HTO")
plot3
dev.off()

# Add ratio calculation between nCount and nFeature. 
L2@meta.data$RNA_ratio <- L2@meta.data$nCount_RNA/L2@meta.data$nFeature_RNA
L2@meta.data$HTO_ratio <- L2@meta.data$nCount_HTO/L2@meta.data$nFeature_HTO

# Removing unwanted cells from the data based on above QC plots
L2 <- subset(L2, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 50)

# Normalize HTOs and subset singlets
L2 <- NormalizeData(L2, assay = "HTO", normalization.method = "CLR")
L2 <- HTODemux(L2, assay = "HTO", positive.quantile = 0.99)
meta <- L2@meta.data
table(meta$HTO_classification.global)

# Group cells based on the max HTO signal
DefaultAssay(L2) <- "HTO"
Idents(L2) <- "HTO_maxID"
pdf(file.path(figures_dir, "QC_RidgePlot_HTO1-8.pdf"),width=10,height=8,paper='special')
RidgePlot(L2, assay = "HTO", features = rownames(L2[["HTO"]])[1:4], ncol = 2)
dev.off()

pdf(file.path(figures_dir, "QC_Scatter_HTO.pdf"),width=10,height=8,paper='special')
plot1 <- FeatureScatter(L2, feature1 = "M-HTO-1", feature2 = "M-HTO-2")
plot2 <- FeatureScatter(L2, feature1 = "M-HTO-3", feature2 = "M-HTO-4")
plot3 <- FeatureScatter(L2, feature1 = "M-HTO-1", feature2 = "M-HTO-3")
plot4 <- FeatureScatter(L2, feature1 = "M-HTO-2", feature2 = "M-HTO-4")
plot1 + plot2 + plot3 + plot4
dev.off()

# subset singlets
Idents(L2) <- "HTO_classification.global"
L2.singlet <- subset(L2, idents = "Singlet")
meta.L2.singlet <- L2.singlet@meta.data

L2.CellBC_WL <- rownames(L2.singlet@meta.data)
L2.CellBC_WL <- sub("-1$", "", L2.CellBC_WL) # Remove the trailing "-1"

writeLines(L2.CellBC_WL, con = "./DATA/BioHPC.files/L2_CellBC_WL.txt")  # this is the file for 03_correctBC_parallel_v3.py

saveRDS(L2.singlet, file = "./DATA/RDS/L2.singlet.rds")

