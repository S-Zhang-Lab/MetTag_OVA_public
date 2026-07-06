
#13.3 cell cycle regression

### Step7.3.2 analysis UTSW02 - cancer cell focused ###

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

Mets <- readRDS("Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")

Mets@meta.data <- Mets@meta.data %>%
  mutate(
    LARRY_class = case_when(
      
      ## ---- om cells ----
      grepl("_Om$", orig.ident) & is.na(om_LARRY_bc_classification) ~ "om_negative",
      grepl("_Om$", orig.ident) & om_LARRY_bc_classification == "singlet"  ~ "om_singlet",
      grepl("_Om$", orig.ident) & om_LARRY_bc_classification == "doublet"  ~ "om_doublet",
      grepl("_Om$", orig.ident) & om_LARRY_bc_classification %in% c("multilet", "multiplet") ~ "om_multiplet",
      
      ## ---- Asc cells ----
      grepl("_Asc$", orig.ident) & is.na(Asc_LARRY_bc_classification) ~ "Asc_negative",
      grepl("_Asc$", orig.ident) & Asc_LARRY_bc_classification == "singlet"  ~ "Asc_singlet",
      grepl("_Asc$", orig.ident) & Asc_LARRY_bc_classification == "doublet"  ~ "Asc_doublet",
      grepl("_Asc$", orig.ident) & Asc_LARRY_bc_classification %in% c("multilet", "multiplet") ~ "Asc_multiplet",
      
      TRUE ~ "unclassified"
    )
  )

Mets$LARRY_class <- factor(
  Mets$LARRY_class,
  levels = c(
    "om_negative", "om_singlet", "om_doublet", "om_multiplet",
    "Asc_negative", "Asc_singlet", "Asc_doublet", "Asc_multiplet",
    "unclassified"
  )
)

table(Mets$LARRY_class, Mets$orig.ident, useNA = "ifany")

DimPlot(
  Mets,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "LARRY_class",
  label = FALSE,
  ncol = 4
)

### REGRESS CELL CYCLE
library(Seurat)

# Save pre-regression UMAP
Mets[["umap_pre_cc"]] <- Mets[["umap"]]

# Load human cell cycle genes
data(cc.genes)

# Convert human gene symbols to mouse style (first letter uppercase)
human_to_mouse <- function(genes) {
  paste0(toupper(substr(genes, 1, 1)), tolower(substr(genes, 2, nchar(genes))))
}

s.genes <- human_to_mouse(cc.genes$s.genes)
g2m.genes <- human_to_mouse(cc.genes$g2m.genes)

# Use SCT assay for integrated dataset
assay_to_use <- "SCT"

# Only keep genes actually present in SCT assay
s.genes <- s.genes[s.genes %in% rownames(Mets[[assay_to_use]])]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(Mets[[assay_to_use]])]

length(s.genes)   # should now be >0
length(g2m.genes) # should now be >0

# Cell cycle scoring
Mets <- CellCycleScoring(
  object = Mets,
  assay = assay_to_use,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# Regress out cell cycle using SCTransform
Mets <- SCTransform(
  object = Mets,
  assay = assay_to_use,
  vars.to.regress = c("S.Score", "G2M.Score"),
  verbose = FALSE
)

# Recompute PCA and UMAP
Mets <- RunPCA(Mets, verbose = FALSE)
Mets <- RunUMAP(Mets, dims = 1:30)

# Compare UMAPs before vs after regression
DimPlot(Mets, reduction = "umap_pre_cc", group.by = "Phase") + ggtitle("Before Cell Cycle Regression")
DimPlot(Mets, reduction = "umap", group.by = "Phase") + ggtitle("After Cell Cycle Regression")

library(patchwork)  # for side-by-side plotting

p1 <- DimPlot(Mets, reduction = "umap_pre_cc", group.by = "Phase") + ggtitle("Before Cell Cycle Regression")
p2 <- DimPlot(Mets, reduction = "umap", group.by = "Phase") + ggtitle("After Cell Cycle Regression")

p1 + p2  # displays them side by side

# Visualize the UMAP after cell cycle regression
DimPlot(Mets, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("UMAP after Cell Cycle Regression (Seurat Clusters)")

DimPlot(Mets, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("UMAP after Cell Cycle Regression (Labeled Clusters)")


# Find markers for all clusters
all_markers <- FindAllMarkers(
  object = Mets,
  only.pos = TRUE,        # keep only genes upregulated in the cluster
  min.pct = 0.25,         # gene must be expressed in at least 25% of cells
  logfc.threshold = 0.25  # minimum log2 fold change
)

# Inspect top markers
head(all_markers)

# Optional: top 10 markers per cluster
top10 <- all_markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 10)

# View top10
top10
# Assuming you already have all_markers from FindAllMarkers()
top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  pull(gene)  # extract gene names

# Plot heatmap
DoHeatmap(Mets, features = top10, size = 3) + 
  ggtitle("Top 10 Marker Genes per Cluster") +
  theme(axis.text.y = element_text(size = 7))

# === Examine IFNa and IFNg genesets ================
# Read hallmark gene set
hallmark_gmt <- read.gmt("mh.all.v2024.1.Mm.symbols.gmt.txt")

# Extract genes for the INTERFERON_ALPHA_RESPONSE set
ifna_genes <- hallmark_gmt %>%
  filter(term == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene)

# View
print(ifna_genes)

# Use your ifna_genes list and intersect with canonical IFN-gamma genes
canonical_ifng_genes <- c(
  # Key transcription factors and signaling
  "Irf1", "Irf8", "Stat1", "Stat2", "Socs1", "Ciita",
  
  # Chemokines
  "Cxcl9", "Cxcl10", "Cxcl11", "Cxcl13",
  
  # Antigen processing & presentation machinery (MHC I & II)
  "Tap1", "Tap2", "Tapbp",  # TAP and associated proteins
  "Psmb8", "Psmb9", "Psme1", "Psme2",  # Immunoproteasome components
  "H2-K1", "H2-D1", "H2-Aa", "H2-Ab1", "H2-Eb1",
  "Cd74", "B2m",
  
  # Guanylate-binding proteins (GBPs)
  "Gbp2b", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Gbp6", "Gbp7",
  
  # Interferon-stimulated genes (ISGs)
  "Isg15", "Ifi30", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oas2",
  "Oas3", "Oasl1", "Rsad2", "Irf7", "Lgals3bp", "Trim21",
  
  # Enzymes and inflammatory mediators
  "Nos2",  # inducible nitric oxide synthase
  "Casp1", "Casp8",  # inflammatory/apoptotic caspases
  "Tnf",  # sometimes induced downstream of IFN-γ
  
  # Immune checkpoint and receptors
  "Cd274",  # PD-L1
  "Ifngr1", "Ifngr2"  # IFN-γ receptor subunits
)


# Overlap
overlap <- intersect(ifna_genes, canonical_ifng_genes)
length(overlap)
print(overlap)

Idents(Mets) <- "Combined.HTO_group"
DimPlot(Mets)
VlnPlot(Mets, features = overlap, assay = "SCT")

pdf(file = "Mets_VlnPlot_IFNa-IFNg_overlap_Asc.VS.Om.pdf", width = 9, height = 18, paper='special')
VlnPlot(Mets, features = overlap, assay = "SCT")
dev.off()

ifna_genes <- tolower(ifna_genes)
canonical_ifng_genes <- tolower(canonical_ifng_genes)
ifna_genes <- ifna_genes[ifna_genes %in% rownames(Mets[["RNA"]])]
canonical_ifng_genes <- canonical_ifng_genes[canonical_ifng_genes %in% rownames(Mets[["RNA"]])]

length(ifna_genes)
length(canonical_ifng_genes)

Mets <- AddModuleScore(
  Mets, 
  features = list(ifna_genes), 
  name = "IFNa_score", 
  assay = "SCT"
)

Mets <- AddModuleScore(
  Mets, 
  features = list(canonical_ifng_genes), 
  name = "IFNG_score", 
  assay = "SCT"
)

VlnPlot(Mets, features = "IFNa_score1", group.by = "Combined.HTO_group")
VlnPlot(Mets, features = "IFNG_score1", group.by = "Combined.HTO_group")
FeaturePlot(Mets, features = c("IFNa_score1", "IFNG_score1"))

# Create a data frame for IFNa
scores_df <- data.frame(
  score = Mets$IFNg_score1,
  group = Mets$Combined.HTO_group
)

# Wilcoxon rank-sum test
wilcox_test <- wilcox.test(score ~ group, data = scores_df)
wilcox_test$p.value

# [1] 7.449182e-93
# 
# t-test
t_test <- t.test(score ~ group, data = scores_df)
t_test$p.value


# Wilcoxon rank-sum test
wilcox_test <- wilcox.test(score ~ group, data = scores_df)
wilcox_test$p.value

# Optional: t-test
t_test <- t.test(score ~ group, data = scores_df)
t_test$p.value

