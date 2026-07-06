# ======================================================
# 1️⃣  Load packages
# ======================================================
library(DESeq2)
library(dplyr)
library(data.table)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)

setwd("/Users/s438978/Desktop/RData/RNA-seq")

# ======================================================
# 2️⃣  Define count files — ULA vs Reg, NTC only
# ======================================================
count_files <- list(
  Reg_NTC_PBS_1 = "Reg_NP_1_counts.txt",
  Reg_NTC_PBS_2 = "Reg_NP_2_counts.txt",
  Reg_NTC_PBS_3 = "Reg_NP_3_counts.txt",
  ULA_NTC_PBS_1 = "ULA_NP_1_counts.txt",
  ULA_NTC_PBS_2 = "ULA_NP_2_counts.txt",
  ULA_NTC_PBS_3 = "ULA_NP_3_counts.txt"
)

# ======================================================
# 3️⃣  Read and combine counts
# ======================================================
count_list <- lapply(count_files, fread, skip = 1)
gene_ids <- count_list[[1]]$Geneid
count_list <- lapply(count_list, function(x) x[[ncol(x)]])
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)

# ======================================================
# 4️⃣  Metadata and DESeq2 setup
# ======================================================
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c(rep("Reg", 3), rep("ULA", 3)))
)
sample_info$condition <- relevel(sample_info$condition, ref = "Reg")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# ======================================================
# 5️⃣  Differential expression (ULA vs Reg)
# ======================================================
resLFC <- lfcShrink(dds, coef = "condition_ULA_vs_Reg", type = "apeglm")
resLFC <- as.data.frame(resLFC[order(resLFC$padj), ])
resLFC$ENSEMBL <- sub("\\..*$", "", rownames(resLFC))

# ======================================================
# 6️⃣  Map Ensembl → gene symbol
# ======================================================
gene_map <- bitr(
  resLFC$ENSEMBL,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Mm.eg.db
)
resLFC$SYMBOL <- gene_map$SYMBOL[match(resLFC$ENSEMBL, gene_map$ENSEMBL)]
resLFC$SYMBOL[is.na(resLFC$SYMBOL)] <- resLFC$ENSEMBL

# ======================================================
# 7️⃣  Volcano plot
# ======================================================
EnhancedVolcano(
  resLFC,
  lab = resLFC$SYMBOL,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "ULA vs Reg (PBS) — NTC only",
  subtitle = "Positive LFC = Up in ULA",
  caption = "padj < 0.05, |LFC| ≥ 1",
  col = c("gray30", "gray70", "skyblue3", "firebrick2")
)

# ======================================================
# 8️⃣  Heatmap of top 50 DEGs
# ======================================================
sig_genes <- rownames(resLFC)[!is.na(resLFC$padj) & resLFC$padj < 0.05]
top50 <- head(sig_genes[order(resLFC$padj[sig_genes])], 50)

vst_mat <- assay(vst(dds, blind = FALSE))
heatmat <- vst_mat[top50, ]
rownames(heatmat) <- resLFC$SYMBOL[match(rownames(heatmat), rownames(resLFC))]

annotation_col <- data.frame(Condition = colData(dds)$condition)
rownames(annotation_col) <- colnames(heatmat)

pheatmap(
  heatmat,
  annotation_col = annotation_col,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  main = "Top 50 DEGs — ULA vs Reg (PBS)"
)

# ======================================================
# 9️⃣  Reactome pathway enrichment — UP & DOWN
# ======================================================

## --- Separate gene sets ---
sig_up <- resLFC %>% filter(padj < 0.05 & log2FoldChange > 0)
sig_down <- resLFC %>% filter(padj < 0.05 & log2FoldChange < 0)

## --- Map to Entrez IDs ---
up_entrez <- bitr(sig_up$ENSEMBL, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
down_entrez <- bitr(sig_down$ENSEMBL, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

## --- Enrich Reactome pathways ---
reactome_up <- enrichPathway(gene = up_entrez$ENTREZID, organism = "mouse", readable = TRUE)
reactome_down <- enrichPathway(gene = down_entrez$ENTREZID, organism = "mouse", readable = TRUE)

## --- Plot top 15 each ---
dotplot(reactome_up, showCategory = 15, title = "Reactome — Up in ULA (vs Reg)")
dotplot(reactome_down, showCategory = 15, title = "Reactome — Up in Reg (vs ULA)")

# ======================================================
# 🔟  Save DE results
# ======================================================
write.csv(resLFC, "DEGs_ULA_vs_Reg_NTC.csv")

## --- Prepare data for bar chart ---
plot_reactome_bar <- function(reactome_result, title_text, color_fill) {
  if (nrow(reactome_result) == 0) return(NULL)
  df <- as.data.frame(reactome_result) %>%
    dplyr::mutate(
      Description = factor(Description, levels = rev(Description)),
      logFDR = -log10(p.adjust)
    ) %>%
    head(15)
  
  ggplot(df, aes(x = Description, y = logFDR, fill = logFDR)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "skyblue1", high = color_fill) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      x = NULL,
      y = expression(-log[10]("FDR")),
      title = title_text
    )
}

# --- Plot ---
plot_reactome_bar(reactome_up, "Reactome Pathways — Up in ULA (vs Reg)", "firebrick3")
plot_reactome_bar(reactome_down, "Reactome Pathways — Up in Reg (vs ULA)", "steelblue4")

# ======================================================

### You can define gene lists manually or download from MSigDB/Seurat
# --- Define gene lists ---
# vst_mat already has unique rownames
# -----------------------------

# -----------------------------
# 1. Define gene sets (mouse)
# -----------------------------
library(dplyr)
library(ggplot2)
library(pheatmap)

# --- 1. Define S and G2/M genes (mouse canonical) ---
S_genes   <- c("Pcna", "Cdk2")  # expand as needed
G2M_genes <- c("Cdc20","Ccnb1","Ccnb2","Cdk1","Top2a","Aurkb","Plk1","Birc5")  # canonical G2/M

cell_cycle_genes <- c(S_genes, G2M_genes)

# --- 2. Prepare VST matrix with gene symbols ---
vst_mat <- assay(vst(dds, blind = FALSE))
gene_symbols <- resLFC$SYMBOL[match(rownames(vst_mat), rownames(resLFC))]

# ensure character type
gene_symbols <- as.character(gene_symbols)

# replace NAs or empty with the original Ensembl IDs
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(vst_mat)[is.na(gene_symbols) | gene_symbols == ""]

# ensure unique and clean rownames
rownames(vst_mat) <- make.unique(gene_symbols)


# Keep only genes in our list
cc_genes_present <- cell_cycle_genes[cell_cycle_genes %in% rownames(vst_mat)]
vst_cc <- vst_mat[cc_genes_present, ]

# --- 3. Scale genes (Z-score by row) ---
vst_cc_scaled <- t(scale(t(vst_cc)))  # scale each gene across samples

# --- 4. Create annotation for samples ---
sample_annotation <- data.frame(
  Condition = colData(dds)$condition
)
rownames(sample_annotation) <- colnames(vst_cc_scaled)

# --- 5. Plot heatmap ---
pheatmap(
  vst_cc_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = sample_annotation,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  main = "S & G2/M Cell Cycle Genes: ULA vs Reg"
)

# === Check "Signaling by TGFB family members" genes and compare ULA vs Reg ===
# === TGFB pathway heatmap: ULA vs Reg ===
# === Cytokine and TGFB receptor expression: ULA vs Reg ===
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)

# --- 1️⃣ Define receptor genes of interest ---
custom_receptors <- c(
  "Tgfbr1", "Tgfbr2", "Tgfbr3",        # TGFB receptors
  "Tnfrsf1a", "Tnfrsf1b",              # TNF receptors
  "Il1r1", "Il1r2", "Il1rap",          # IL1 receptors
  "Il10ra", "Il10rb"                   # IL10 receptors
)

# --- 2️⃣ Choose your VST matrix ---
# You probably have one of these available; adjust the name if needed:
# vst_scaled <- t(scale(t(vst_sym)))   # if you only have vst_sym
# or if vst_scaled already exists, just use that

# confirm which matrix exists
if (exists("vst_sym")) {
  vst_mat <- vst_sym
} else if (exists("vst_scaled")) {
  vst_mat <- vst_scaled
} else if (exists("vst_mat")) {
  vst_mat <- vst_mat
} else {
  stop("No VST matrix found (expected one of vst_sym, vst_scaled, or vst_mat).")
}

# --- 3️⃣ Keep only receptors present in the dataset ---
receptors_present <- intersect(custom_receptors, rownames(vst_mat))
if (length(receptors_present) < 2)
  stop("Too few receptor genes found in VST matrix. Check identifiers!")

# --- 4️⃣ Subset VST matrix and scale rows (Z-score) ---
vst_receptors <- vst_mat[receptors_present, , drop = FALSE]
vst_receptors_scaled <- t(scale(t(vst_receptors)))

# --- 5️⃣ Create sample annotation ---
annotation_col <- data.frame(Condition = colData(dds)$condition)
rownames(annotation_col) <- colnames(vst_receptors_scaled)

# --- 6️⃣ Plot heatmap ---
pheatmap(vst_receptors_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         main = "TGFB and Cytokine Receptor Expression (ULA vs Reg)",
         fontsize_row = 10,
         show_rownames = TRUE)

# --- 7️⃣ Compute summary stats per receptor ---
expr_long <- vst_receptors_scaled %>%
  as.data.frame() %>%
  rownames_to_column(var = "SYMBOL") %>%
  pivot_longer(cols = -SYMBOL, names_to = "Sample", values_to = "VST") %>%
  left_join(
    data.frame(Sample = colnames(vst_receptors_scaled),
               Condition = as.character(colData(dds)$condition),
               stringsAsFactors = FALSE),
    by = "Sample"
  )

gene_stats <- expr_long %>%
  group_by(SYMBOL) %>%
  summarize(
    mean_Reg = mean(VST[Condition == "Reg"], na.rm = TRUE),
    mean_ULA = mean(VST[Condition == "ULA"], na.rm = TRUE),
    med_Reg = median(VST[Condition == "Reg"], na.rm = TRUE),
    med_ULA = median(VST[Condition == "ULA"], na.rm = TRUE),
    wilcox_p = tryCatch(wilcox.test(VST ~ Condition, data = cur_data_all())$p.value,
                        error = function(e) NA_real_),
    n_Reg = sum(Condition == "Reg"),
    n_ULA = sum(Condition == "ULA"),
    .groups = "drop"
  ) %>%
  mutate(wilcox_p_adj = p.adjust(wilcox_p, method = "BH"))

# --- 8️⃣ View per-gene stats ---
print(gene_stats)

