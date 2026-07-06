# -----------------------------
# 0. Load libraries
# -----------------------------
library(data.table)
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(dplyr)
library(RColorBrewer)

# -----------------------------
# 1. Set working directory
# -----------------------------
setwd("/Users/s438978/Desktop/RData/RNA-seq")

# -----------------------------
# 2. Define count files
# -----------------------------
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

# -----------------------------
# 3. Read and combine counts
# -----------------------------
count_list <- lapply(count_files, function(f) {
  dt <- fread(f, skip = 1)
  if (ncol(dt) < 2) stop(paste("File has too few columns:", f))
  dt
})

# Extract gene IDs (from first file)
gene_ids <- count_list[[1]]$Geneid

# Extract counts (last column of each file)
count_list <- lapply(count_list, function(x) x[[ncol(x)]])
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)

# -----------------------------
# 4. Sample metadata
# -----------------------------
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  plate = rep(c("Reg", "ULA"), each = 12),
  genotype = rep(c("NTC", "KO"), each = 6, times = 2),
  treatment = rep(c("PBS", "IFNg"), each = 3, times = 4)
)

# Combine into a single condition factor (optional)
sample_info$condition <- factor(paste(sample_info$plate, sample_info$genotype, sample_info$treatment, sep = "_"))

# -----------------------------
# 5. Create DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

# -----------------------------
# 6. Variance-stabilizing transformation
# -----------------------------
vst_mat <- assay(vst(dds, blind = FALSE))

# -----------------------------
# 7. Map gene symbols
# -----------------------------
ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

gene_map <- getBM(
  attributes=c("ensembl_gene_id","mgi_symbol"),
  filters="ensembl_gene_id",
  values=rownames(dds),
  mart=ensembl
)
gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]

gene_symbols <- gene_map$mgi_symbol[match(rownames(dds), gene_map$ensembl_gene_id)]
gene_symbols[is.na(gene_symbols) | gene_symbols==""] <- rownames(dds)
names(gene_symbols) <- rownames(dds)

# -----------------------------
# 8. Define TNF pathway genes
# -----------------------------
canonical_tnf_genes <- c(
  "Tnf", "Tnfrsf1a", "Tnfrsf1b", "Traf2", "Traf5", "Ripk1", "Tab2", "Ikbkg",
  "Ikbkb", "Nfkb1", "Nfkb2", "Rel", "Relb", "Casp8", "Casp3", "Cflar",
  "Birc2", "Xiap", "Smpd2", "Smpd3", "Ubc"
)

tnf_genes_present <- gene_symbols[names(gene_symbols) %in% rownames(vst_mat) &
                                    gene_symbols %in% canonical_tnf_genes]

tnf_vst_mat <- vst_mat[names(tnf_genes_present), ]
rownames(tnf_vst_mat) <- tnf_genes_present

# -----------------------------
# 9. Row annotation (anti vs pro-apoptotic)
# -----------------------------
anti_apoptotic <- c("Xiap", "Birc2", "Cflar", "Ikbkb", "Ikbkg", "Tab2", "Traf5", "Rel")
pro_apoptotic  <- c("Tnfrsf1a", "Casp8", "Casp3", "Smpd2", "Smpd3")

row_annotation <- data.frame(
  Function = ifelse(rownames(tnf_vst_mat) %in% anti_apoptotic, "Anti-apoptotic", "Pro-apoptotic")
)
rownames(row_annotation) <- rownames(tnf_vst_mat)

# -----------------------------
# 10. Column annotation and order by treatment
# -----------------------------
col_annotation <- data.frame(
  Plate = sample_info$plate,
  Genotype = sample_info$genotype,
  Treatment = sample_info$treatment
)
rownames(col_annotation) <- rownames(sample_info)

# Order columns: PBS first, then IFNg
col_order <- order(col_annotation$Treatment, col_annotation$Plate, col_annotation$Genotype)
tnf_vst_mat <- tnf_vst_mat[, col_order]
col_annotation <- col_annotation[col_order, ]

# -----------------------------
# 11. Plot TNF pathway heatmap
# -----------------------------
pheatmap(
  tnf_vst_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # preserve manual ordering by treatment
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "TNF Pathway Genes Across All Groups (PBS → IFNg)"
)

# -----------------------------
# 1. Subset count matrix and sample info
# -----------------------------
groups_of_interest <- c("ULA_NTC_PBS", "ULA_NTC_IFNg")
sample_mask <- grepl("ULA_NTC_PBS|ULA_NTC_IFNg", colnames(count_matrix))
count_sub <- count_matrix[, sample_mask]
sample_info_sub <- sample_info[sample_mask, , drop = FALSE]

# Add a clean 'condition' column
sample_info_sub$condition <- factor(ifelse(grepl("PBS", rownames(sample_info_sub)), "PBS", "IFNg"),
                                    levels = c("PBS", "IFNg"))

# -----------------------------
# 2. Create DESeq2 object
# -----------------------------
dds_sub <- DESeqDataSetFromMatrix(
  countData = count_sub,
  colData = sample_info_sub,
  design = ~ condition
)

# Pre-filter low-count genes
dds_sub <- dds_sub[rowSums(counts(dds_sub)) >= 10, ]

# Run DESeq2
dds_sub <- DESeq(dds_sub)

# -----------------------------
# 3. Extract results: IFNg vs PBS
# -----------------------------
res_sub <- results(dds_sub, contrast = c("condition", "IFNg", "PBS"))

# Add gene symbols
ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
gene_map_sub <- getBM(
  attributes=c("ensembl_gene_id","mgi_symbol"),
  filters="ensembl_gene_id",
  values=rownames(res_sub),
  mart=ensembl
)
gene_map_sub <- gene_map_sub[!duplicated(gene_map_sub$ensembl_gene_id), ]
res_sub$gene_symbol <- gene_map_sub$mgi_symbol[match(rownames(res_sub), gene_map_sub$ensembl_gene_id)]
res_sub$gene_symbol[is.na(res_sub$gene_symbol)] <- rownames(res_sub)

# -----------------------------
# 4. Save DEGs table
# -----------------------------
res_sub_df <- as.data.frame(res_sub)
res_sub_df <- res_sub_df[order(res_sub_df$padj), ]
write.csv(res_sub_df, "DEGs_ULANC_IFNg_vs_PBS.csv", row.names = TRUE)

# -----------------------------
# 5. Volcano plot
# -----------------------------
library(EnhancedVolcano)
EnhancedVolcano(
  res_sub_df,
  lab = res_sub_df$gene_symbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "ULANC IFNg vs PBS",
  subtitle = "NTC cells, ULA plate",
  drawConnectors = TRUE,
  widthConnectors = 0.5
)

# -----------------------------
# 6. Heatmap of top 50 DEGs
# -----------------------------
library(pheatmap)
vst_sub <- assay(vst(dds_sub, blind=FALSE))
top_genes <- rownames(res_sub_df)[1:50]  # top 50 by adjusted p-value
heatmat <- vst_sub[top_genes, ]
rownames(heatmat) <- res_sub_df$gene_symbol[match(rownames(heatmat), rownames(res_sub_df))]

# Correct annotation for pheatmap
annotation_col <- data.frame(Treatment = sample_info_sub$condition)
rownames(annotation_col) <- colnames(heatmat)  # must match columns

pheatmap(
  heatmat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Top 50 DEGs: ULA NTC IFNg vs PBS"
)

# =====================================================
# 🧬 RNA-seq DESeq2: ULA and Reg NTC (IFNg vs PBS)
# =====================================================

# -----------------------------
# 0. Load libraries
# -----------------------------
library(data.table)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(EnhancedVolcano)
library(patchwork)   # for combining plots

# -----------------------------
# 1. Set working directory
# -----------------------------
setwd("/Users/s438978/Desktop/RData/RNA-seq")

# -----------------------------
# 2. Define count files
# -----------------------------
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

# -----------------------------
# 3. Read and combine counts
# -----------------------------
count_list <- lapply(count_files, function(f) {
  dt <- fread(f, skip = 1)
  if (ncol(dt) < 2) stop(paste("File has too few columns:", f))
  dt
})

# Extract gene IDs
gene_ids <- count_list[[1]]$Geneid
count_list <- lapply(count_list, function(x) x[[ncol(x)]])
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)

# -----------------------------
# 4. Sample metadata
# -----------------------------
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  plate = rep(c("Reg", "ULA"), each = 12),
  genotype = rep(c("NTC", "KO"), each = 6, times = 2),
  treatment = rep(c("PBS", "IFNg"), each = 3, times = 4)
)
sample_info$condition <- factor(paste(sample_info$plate, sample_info$genotype, sample_info$treatment, sep = "_"))

# -----------------------------
# 5. Create DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# =====================================================
# 🔹 Function: Run DESeq2 contrast and volcano plot
# =====================================================
run_volcano <- function(count_matrix, sample_info, pattern, output_csv, title, subtitle) {
  # Subset for desired samples
  mask <- grepl(pattern, colnames(count_matrix))
  count_sub <- count_matrix[, mask]
  sample_info_sub <- sample_info[mask, , drop = FALSE]
  
  sample_info_sub$condition <- factor(
    ifelse(grepl("PBS", rownames(sample_info_sub)), "PBS", "IFNg"),
    levels = c("PBS", "IFNg")
  )
  
  # DESeq2
  dds_sub <- DESeqDataSetFromMatrix(count_sub, sample_info_sub, design = ~ condition)
  dds_sub <- dds_sub[rowSums(counts(dds_sub)) >= 10, ]
  dds_sub <- DESeq(dds_sub)
  res_sub <- results(dds_sub, contrast = c("condition", "IFNg", "PBS"))
  
  # Map gene symbols
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(res_sub),
    mart = ensembl
  )
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  res_sub$gene_symbol <- gene_map$mgi_symbol[match(rownames(res_sub), gene_map$ensembl_gene_id)]
  res_sub$gene_symbol[is.na(res_sub$gene_symbol)] <- rownames(res_sub)
  
  # Save DEGs
  res_df <- as.data.frame(res_sub)
  res_df <- res_df[order(res_df$padj), ]
  write.csv(res_df, output_csv, row.names = TRUE)
  
  # Volcano plot
  p <- EnhancedVolcano(
    res_df,
    lab = res_df$gene_symbol,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    title = title,
    subtitle = subtitle,
    drawConnectors = TRUE,
    widthConnectors = 0.5
  )
  
  return(p)
}

# =====================================================
# 🔹 Run comparisons
# =====================================================

# (1) ULA NTC IFNg vs PBS
p_ula <- run_volcano(
  count_matrix,
  sample_info,
  pattern = "ULA_NTC_PBS|ULA_NTC_IFNg",
  output_csv = "DEGs_ULANTC_PBS_vs_IFNg.csv",
  title = "ULA NTC: IFNg vs PBS",
  subtitle = "NTC cells, ULA plate"
)

# (2) Reg NTC IFNg vs PBS
p_reg <- run_volcano(
  count_matrix,
  sample_info,
  pattern = "Reg_NTC_PBS|Reg_NTC_IFNg",
  output_csv = "DEGs_RegNTC_PBS_vs_IFNg.csv",
  title = "Reg NTC: IFNg vs PBS",
  subtitle = "NTC cells, Reg plate"
)

# =====================================================
# 🔹 Combine volcano plots side-by-side
# =====================================================
p_combined <- p_reg + p_ula + plot_layout(ncol = 2)
p_combined



