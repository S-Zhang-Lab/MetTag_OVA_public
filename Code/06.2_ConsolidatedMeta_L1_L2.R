# === Consolidate L1-L2 data =====

# Load Packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# load L1-L2 transcriptome and HTO data from cellranger output
base_dir <- "C:/Users/s208205/Downloads/R_projects/MetTag_BC_OVA"
figures_dir <- file.path(base_dir, "Combined_Figures")

L1 <- readRDS("./DATA/RDS/L1.singlet_post_BC_assignement.rds")
L2 <- readRDS("./DATA/RDS/L2.singlet_post_BC_assignement.rds")

# Rename col names in meta.data
colnames(L1@meta.data) <- paste0("om_", colnames(L1@meta.data))
colnames(L2@meta.data) <- paste0("Asc_", colnames(L2@meta.data))

colnames(L1@meta.data)
colnames(L2@meta.data)

## recode HTO groups ##
L1@meta.data$om_HTO_group <- recode(L1@meta.data$om_hash.ID, 
                                                "M-HTO-5" = "4_NaiveOm", 
                                                "M-HTO-3" = "2_OmMet", 
                                                "M-HTO-4" = "2_OmMet")

 
L2@meta.data$Asc_HTO_group <- recode(L2@meta.data$Asc_hash.ID, 
                                    "M-HTO-5" = "4_NaiveOm", 
                                    "M-HTO-3" = "2_OmMet", 
                                    "M-HTO-4" = "2_OmMet"
                                    )


L2@meta.data$Asc_HTO_group <- recode(L2@meta.data$Asc_hash.ID, 
                                      "M-HTO-1" = "3_Lavage", 
                                      "M-HTO-2" = "3_Lavage", 
                                      "M-HTO-3" = "3_Lavage", 
                                      "M-HTO-4" = "3_Lavage", 
                                      "M-HTO-5" = "1_AscMet", 
                                      "M-HTO-6" = "1_AscMet", 
                                      "M-HTO-7" = "1_AscMet", 
                                      "M-HTO-8" = "1_AscMet"
                                     )


meta.L1 <- L1@meta.data
meta.L2 <- L2@meta.data

# ======   merge the two Seurat objects ==========
Ova_merged <- merge(x = L1, y = L2, add.cell.ids = c("Om", "Asc"), project = "Ova_Merged_Data")
meta <- Ova_merged@meta.data

# consolidate meta further by combining HTO from L1 and L2 as new "Combined.HTO_group") 
colnames(meta)
df <- meta
df <- df %>%
  mutate(Combined.HTO_group = ifelse(!is.na(om_HTO_group) & om_HTO_group != "", 
                                     paste0("om:", om_HTO_group), 
                                     ifelse(!is.na(Asc_HTO_group) & Asc_HTO_group != "", 
                                            paste0("asc:", Asc_HTO_group), 
                                            NA)))
table(df$Combined.HTO_group)
Ova_merged@meta.data <- df

####  Consolidate LARRY BC annotation and check any common LARRY BC in both om and asc samples
meta <- Ova_merged@meta.data
colnames(meta)
df <- meta
df <- df %>%
  mutate(Combined.LARRY_group = ifelse(!is.na(om_LARRY_first_bc) & om_LARRY_first_bc != "", 
                                       paste0(om_LARRY_first_bc, "-om"), 
                                       ifelse(!is.na(Asc_LARRY_first_bc) & Asc_LARRY_first_bc != "", 
                                              paste0(Asc_LARRY_first_bc,"-asc"), 
                                              NA)))
table(df$Combined.LARRY_group)

df <- df %>%
  mutate(orig.ident = ifelse(!is.na(om_orig.ident.x) & om_orig.ident.x != "", 
                             om_orig.ident.x, 
                            ifelse(!is.na(Asc_orig.ident.x) & Asc_orig.ident.x != "", 
                            Asc_orig.ident.x, 
                            NA)))

table(df$orig.ident)

# write back
Ova_merged@meta.data <- df 


###  ===== check merged obj structure  =======
DefaultAssay(Ova_merged) <- "RNA"

# we should see "2 layers present: counts.L1_Om, counts.L2_Asc
# 4 other assays present: CITE, HTO, LARRY, Lib_ID
Ova_merged 

### copy obj to for SCT based integration (later)
Ova_merged_SCT <- Ova_merged

## =====  run a standard scRNA-seq analysis before integration ====== ##
# Note: since the data is split into layers, 
# normalization and variable feature identification is performed for each batch independently 
# (a consensus set of variable features is automatically identified)

Ova_merged <- NormalizeData(Ova_merged)
Ova_merged <- FindVariableFeatures(Ova_merged)
Ova_merged <- ScaleData(Ova_merged)
Ova_merged <- RunPCA(Ova_merged)

# now visualize the results of a standard analysis without integration
set.seed(1)
Ova_merged <- FindNeighbors(Ova_merged, dims = 1:30, reduction = "pca")
Ova_merged <- FindClusters(Ova_merged, resolution = 1.5, cluster.name = "unintegrated_clusters")
Ova_merged <- RunUMAP(Ova_merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# visualize by batch and cell type annotation
pdf(file.path(figures_dir, "UMAP_un_integrated.pdf"),width=16,height=5,paper='special')
DimPlot(Ova_merged, reduction = "umap.unintegrated", group.by = c("orig.ident", "om_HTO_group", "Asc_HTO_group"))
dev.off()

### ======  perform integration using various methods  ==== ####

Ova_merged <- IntegrateLayers(
  object = Ova_merged, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE
)

Ova_merged <- IntegrateLayers(
  object = Ova_merged, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)

Ova_merged <- IntegrateLayers(
  object = Ova_merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

Ova_merged <- FindNeighbors(Ova_merged, reduction = "integrated.cca", dims = 1:30)
Ova_merged <- FindClusters(Ova_merged, resolution = 1.5, cluster.name = "cca_clusters")

Ova_merged <- FindNeighbors(Ova_merged, reduction = "integrated.rpca", dims = 1:30)
Ova_merged <- FindClusters(Ova_merged, resolution = 1.5, cluster.name = "rpca_clusters")

Ova_merged <- FindNeighbors(Ova_merged, reduction = "harmony", dims = 1:30)
Ova_merged <- FindClusters(Ova_merged, resolution = 1.5, cluster.name = "harmony_clusters")

# rerun UMAP using integrated data
set.seed(1)
Ova_merged <- RunUMAP(Ova_merged, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
  Ova_merged,
  reduction = "umap.cca",
  group.by = c("orig.ident", "cca_clusters"),
  combine = FALSE, label.size = 2
)

Ova_merged <- RunUMAP(Ova_merged, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p2 <- DimPlot(
  Ova_merged,
  reduction = "umap.rpca",
  group.by = c("orig.ident", "rpca_clusters"),
  combine = FALSE, label.size = 2
)

Ova_merged <- RunUMAP(Ova_merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p3 <- DimPlot(
  Ova_merged,
  reduction = "umap.harmony",
  group.by = c("orig.ident", "harmony_clusters"),
  combine = FALSE, label.size = 2
)

pdf(file.path(figures_dir, "UMAP_post_integration.pdf"),width=19,height=10,paper='special')
wrap_plots(c(p1, p2, p3), ncol = 3, byrow = F)
dev.off()

# rejoin integrated layers (L1 and L2) 
Ova_merged <- JoinLayers(Ova_merged)
Ova_merged

pdf(file.path(figures_dir, "Cca_UMAP_post_integration.pdf"),width=12,height=5,paper='special')
DimPlot(
  Ova_merged,
  reduction = "umap.cca",
  group.by = c("Combined.HTO_group", "cca_clusters"),
  combine = TRUE, label.size = 2
)
dev.off()

# write RDS for downstream analysis. This object has both LARRY BC and LibID annotated 
# can be read later using: Ova_merged <- readRDS("./DATA/RDS/MetTag.ova_L1L2_merged_integrated.rds")
saveRDS(Ova_merged, file = "./DATA/RDS/MetTag.ova_L1L2_merged_integrated.rds")


### =======  alternative integration after SCT transform ======== ####
options(future.globals.maxSize = 3e+09)
Ova_merged_SCT <- SCTransform(Ova_merged_SCT)
Ova_merged_SCT <- RunPCA(Ova_merged_SCT, npcs = 30, verbose = F)
Ova_merged_SCT <- IntegrateLayers(
  object = Ova_merged_SCT,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
Ova_merged_SCT <- FindNeighbors(Ova_merged_SCT, dims = 1:30, reduction = "integrated.dr")
Ova_merged_SCT <- FindClusters(Ova_merged_SCT, resolution = 2)
Ova_merged_SCT <- RunUMAP(Ova_merged_SCT, dims = 1:30, reduction = "integrated.dr")

Ova_merged_SCT

pdf(file.path(figures_dir, "SCT_UMAP_post_integration.pdf"),width=12,height=5,paper='special')
DimPlot(
  Ova_merged_SCT,
  combine = TRUE, label.size = 2, 
  group.by = c("Combined.HTO_group","seurat_clusters")
)
dev.off()

saveRDS(Ova_merged_SCT, file = "./DATA/RDS/MetTag.ova_L1L2_merged_Post_SCT_integrated.rds")

# Read RDS for further analysis
Ova_merged_SCT <- readRDS("./DATA/RDS/MetTag.ova_L1L2_merged_Post_SCT_integrated.rds")
meta.Ova.merged_SCT <- Ova_merged_SCT@meta.data

Ova_merged_SCT_EA <- readRDS("./DATA/EA_SCT/Ova_merged_SCT_EA.rds")
meta.Ova.merged_SCT_EA <- Ova_merged_SCT_EA@meta.data

commmon_cells <- intersect(row.names(meta.Ova.merged_SCT), row.names(meta.Ova.merged_SCT_EA))

colnames(meta.Ova.merged_SCT_EA) 
colnames(meta.Ova.merged_SCT) 

# Add "EA_" prefix to the first table's column names
colnames(meta.Ova.merged_SCT_EA) <- paste("EA", colnames(meta.Ova.merged_SCT_EA), sep = "_")

# Add "New_" prefix to the second table's column names
colnames(meta.Ova.merged_SCT) <- paste("New", colnames(meta.Ova.merged_SCT), sep = "_")

dim(meta.Ova.merged_SCT_EA)
dim(meta.Ova.merged_SCT)

# Convert row names to a new column for merging
meta.Ova.merged_SCT_EA$Cell_names <- rownames(meta.Ova.merged_SCT_EA)
meta.Ova.merged_SCT$Cell_names <- rownames(meta.Ova.merged_SCT)

# Merge data frames by new 'Cell_names' column, keeping all rows from meta.Ova.merged_SCT_EA
merged_data <- merge(meta.Ova.merged_SCT_EA, meta.Ova.merged_SCT, by = "Cell_names", all.x = TRUE)

# View the dimensions and the result
dim(merged_data)
head(merged_data)

row.names(merged_data) <- merged_data$Cell_names
head(merged_data)

# write back to EASCT object.
Ova_merged_SCT_EA@meta.data <- merged_data
saveRDS(Ova_merged_SCT, file = "./DATA/RDS/MetTag.Ova_merged_SCT_with_updated_meta.rds")
