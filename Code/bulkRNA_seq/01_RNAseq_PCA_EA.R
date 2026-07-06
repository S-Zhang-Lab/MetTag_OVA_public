# --------------------------------------------------
# Bulk RNA-seq DESeq2 Analysis — Cleaned & Corrected
# --------------------------------------------------

# 1. Load libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(apeglm)
})

setwd("/Users/s438978/Desktop/RData/RNA-seq")

# 2. List of count files
count_files <- list(
  Reg_NTC_PBS_1   = "Reg_NP_1_counts.txt",
  Reg_NTC_PBS_2   = "Reg_NP_2_counts.txt",
  Reg_NTC_PBS_3   = "Reg_NP_3_counts.txt",
  Reg_NTC_IFNg_1  = "Reg_NI_1_counts.txt",
  Reg_NTC_IFNg_2  = "Reg_NI_2_counts.txt",
  Reg_NTC_IFNg_3  = "Reg_NI_3_counts.txt",
  Reg_KO_PBS_1    = "Reg_C2P_1_counts.txt",
  Reg_KO_PBS_2    = "Reg_C2P_2_counts.txt",
  Reg_KO_PBS_3    = "Reg_C2P_3_counts.txt",
  Reg_KO_IFNg_1   = "Reg_C2I_1_counts.txt",
  Reg_KO_IFNg_2   = "Reg_C2I_2_counts.txt",
  Reg_KO_IFNg_3   = "Reg_C2I_3_counts.txt",
  ULA_NTC_PBS_1   = "ULA_NP_1_counts.txt",
  ULA_NTC_PBS_2   = "ULA_NP_2_counts.txt",
  ULA_NTC_PBS_3   = "ULA_NP_3_counts.txt",
  ULA_NTC_IFNg_1  = "ULA_NI_1_counts.txt",
  ULA_NTC_IFNg_2  = "ULA_NI_2_counts.txt",
  ULA_NTC_IFNg_3  = "ULA_NI_3_counts.txt",
  ULA_KO_PBS_1    = "ULA_C2P_1_counts.txt",
  ULA_KO_PBS_2    = "ULA_C2P_2_counts.txt",
  ULA_KO_PBS_3    = "ULA_C2P_3_counts.txt",
  ULA_KO_IFNg_1   = "ULA_C2I_1_counts.txt",
  ULA_KO_IFNg_2   = "ULA_C2I_2_counts.txt",
  ULA_KO_IFNg_3   = "ULA_C2I_3_counts.txt"
)

# 3. Read and combine counts
count_list <- lapply(count_files, function(f) {
  dt <- fread(f, skip = 1)
  if (ncol(dt) < 2) stop(paste("File has too few columns:", f))
  dt
})

gene_ids <- count_list[[1]]$Geneid
count_list_counts <- lapply(count_list, function(x) x[[ncol(x)]])
count_matrix <- do.call(cbind, count_list_counts)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)

# 4. Metadata
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c(
    rep("Reg_NTC_PBS", 3),
    rep("Reg_NTC_IFNg", 3),
    rep("Reg_KO_PBS", 3),
    rep("Reg_KO_IFNg", 3),
    rep("ULA_NTC_PBS", 3),
    rep("ULA_NTC_IFNg", 3),
    rep("ULA_KO_PBS", 3),
    rep("ULA_KO_IFNg", 3)
  ))
)

# 5. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq normalization and modeling
dds <- DESeq(dds)

# 6. Variance-stabilizing transformation (for QC/PCA)
vst_mat <- vst(dds, blind = TRUE)

# ----------------------------------------------------------
# 7. QC plots (Cook’s distance, dispersion, MA) - FIXED
# ----------------------------------------------------------
pdf("QC_combined_plots_all.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(5, 5, 3, 2))

# 1️⃣ Cook’s Distance
cooks_vals <- assays(dds)[["cooks"]]
if (all(is.na(cooks_vals)) || all(cooks_vals == 0)) {
  plot.new()
  title("Cook’s Distance (no nonzero values)")
} else {
  boxplot(
    log10(cooks_vals + 1e-6),  # add small offset to avoid log(0)
    range = 0,
    las = 2,
    ylab = "log10(Cook’s distance)",
    main = "Cook’s Distance"
  )
}

# 2️⃣ Dispersion Estimates
tryCatch({
  plotDispEsts(dds, main = "Dispersion Estimates")
}, error = function(e) {
  plot.new()
  title("Dispersion plot failed")
  message("Dispersion plot error: ", e$message)
})

# 3️⃣ MA Plot: ULA_NTC_PBS vs Reg_NTC_PBS
res <- results(dds, contrast = c("condition", "ULA_NTC_PBS", "Reg_NTC_PBS"))
resLFC <- lfcShrink(dds, contrast = c("condition", "ULA_NTC_PBS", "Reg_NTC_PBS"), type = "normal")

if (sum(is.finite(resLFC$log2FoldChange)) > 0) {
  plotMA(resLFC, ylim = c(-2, 2), main = "MA Plot: ULA_NTC_PBS vs Reg_NTC_PBS")
} else {
  plot.new()
  title("MA Plot (no finite log2FC values)")
}

dev.off()

# 8. PCA
pcaData <- plotPCA(vst_mat, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = rownames(pcaData)), max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA of 8 Groups (ULA and Reg, KO/NTC, PBS/IFNg)") +
  theme(legend.position = "right")

# Save DE results
write.csv(as.data.frame(resLFC), "DESeq2_results_ULA_NTC_PBS_vs_Reg_NTC_PBS.csv")

# --- 9. Cell cycle heatmap across ULA groups---
library(DESeq2)
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(pheatmap)
# --- Subset VST for ULA samples ---
ula_samples <- grep("^ULA", colnames(vst_mat), value = TRUE)
vst_ula <- assay(vst(dds, blind = TRUE))[, ula_samples]  # <- IMPORTANT: use assay()

# --- Define cell cycle genes ---
G1_genes <- c("Cdk4","Cdk6","Ccnd1","Ccnd2","Ccne1","Ccne2","Rb1","E2f1")
S_genes  <- c("Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7","Pcna","Cdc45")
G2M_genes <- c("Cdk1","Ccnb1","Ccnb2","Cdc20","Cdc25c","Bub1","Plk1","Aurkb")
ordered_genes <- c(G1_genes, S_genes, G2M_genes)

# --- Map symbols to Ensembl IDs ---
symbol2ensembl <- mapIds(org.Mm.eg.db,
                         keys = ordered_genes,
                         column = "ENSEMBL",
                         keytype = "SYMBOL",
                         multiVals = "first")

# --- Keep only genes present ---
genes_present <- symbol2ensembl[symbol2ensembl %in% rownames(vst_ula)]
heatmap_mat <- vst_ula[genes_present, , drop = FALSE]  # numeric matrix

# --- Replace rownames with symbols ---
rownames(heatmap_mat) <- names(genes_present)

# --- Column annotation ---
anno_df <- data.frame(
  Genotype = ifelse(grepl("_KO_", colnames(heatmap_mat)), "KO", "NTC"),
  Treatment = ifelse(grepl("_IFNg_", colnames(heatmap_mat)), "IFNg", "PBS")
)
rownames(anno_df) <- colnames(heatmap_mat)

# --- Order samples ---
sample_order <- anno_df %>%
  arrange(Genotype, Treatment) %>%
  rownames()
heatmap_mat <- heatmap_mat[, sample_order]
anno_df <- anno_df[sample_order, ]

# --- Row annotation (phase) ---
gene_phase <- data.frame(
  Phase = factor(
    ifelse(names(genes_present) %in% G1_genes, "G1",
           ifelse(names(genes_present) %in% S_genes, "S", "G2/M")),
    levels = c("G1","S","G2/M")
  )
)
rownames(gene_phase) <- rownames(heatmap_mat)

# --- Colors ---
ann_colors <- list(
  Genotype = c(NTC = "#1b9e77", KO = "#d95f02"),
  Treatment = c(PBS = "#7570b3", IFNg = "#e7298a")
)
phase_colors <- list(Phase = c(G1 = "#66c2a5", S = "#fc8d62", `G2/M` = "#8da0cb"))


# --- Plot ---
pheatmap(
  heatmap_mat,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = anno_df,
  annotation_row = gene_phase,
  annotation_colors = c(ann_colors, phase_colors),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 7,
  angle_col = 45,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  border_color = NA,
  main = "Cell Cycle Gene Expression (VST) in ULA Samples"
)

# --- 1️⃣ Subset VST matrix to ULA_NTC samples only ---
# LFC shrinkage for IFNg vs PBS in ULA NTC samples
resLFC <- lfcShrink(
  dds,
  contrast = c("condition", "ULA_NTC_IFNg", "ULA_NTC_PBS"),
  type = "normal"   # or 'apeglm' if using a coef
)

# Define cell cycle genes
ordered_genes <- c(G1_genes, S_genes, G2M_genes)

# Map to Ensembl IDs
# --- Clean Ensembl IDs in resLFC ---
rownames(resLFC) <- sub("\\..*", "", rownames(resLFC))  # removes ".1", ".2", etc. if present

# --- Map SYMBOL → ENSEMBL ---
symbol2ensembl <- mapIds(
  org.Mm.eg.db,
  keys = c(G1_genes, S_genes, G2M_genes),
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

# --- Keep only genes present in resLFC ---
present <- names(symbol2ensembl)[symbol2ensembl %in% rownames(resLFC)]

# --- Extract log2FC values ---
lfc_vals <- resLFC[symbol2ensembl[present], "log2FoldChange"]
names(lfc_vals) <- present
lfc_vals <- lfc_vals[!is.na(lfc_vals)]

# --- Make matrix ---
lfc_mat <- matrix(lfc_vals, ncol = 1)
rownames(lfc_mat) <- names(lfc_vals)
colnames(lfc_mat) <- "IFNg_vs_PBS"

# --- Sanity check ---
dim(lfc_mat)
head(lfc_mat)


# Row annotation (cell cycle phase)
gene_phase <- data.frame(
  Phase = factor(
    ifelse(rownames(lfc_mat) %in% G1_genes, "G1",
           ifelse(rownames(lfc_mat) %in% S_genes, "S", "G2/M")),
    levels = c("G1","S","G2/M")
  )
)
rownames(gene_phase) <- rownames(lfc_mat)

# Colors
phase_colors <- list(Phase = c(G1 = "#66c2a5", S = "#fc8d62", `G2/M` = "#8da0cb"))

# Plot
pheatmap(
  lfc_mat,
  scale = "none",           # do not row-normalize for LFC
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = gene_phase,
  annotation_colors = phase_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  border_color = NA,
  main = "Cell Cycle Gene Log2 Fold Change - ULA NTC IFNg vs PBS"
)


