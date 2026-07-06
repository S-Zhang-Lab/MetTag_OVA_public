# =========================================
# 0️⃣ Load libraries
# =========================================
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ReactomePA)
library(msigdbr)
library(patchwork)
library(forcats)
library(biomaRt)
library(pheatmap)

# =========================================
# 1️⃣ Define count files
# =========================================

setwd("/Users/s438978/Desktop/RData/RNA-seq")

count_files <- list(
Reg_NTC_PBS_1   = "Reg_NP_1_counts.txt",
Reg_NTC_PBS_2   = "Reg_NP_2_counts.txt",
Reg_NTC_PBS_3   = "Reg_NP_3_counts.txt",
Reg_NTC_IFNg_1  = "Reg_NI_1_counts.txt",
Reg_NTC_IFNg_2  = "Reg_NI_2_counts.txt",
Reg_NTC_IFNg_3  = "Reg_NI_3_counts.txt",
ULA_NTC_PBS_1   = "ULA_NP_1_counts.txt",
ULA_NTC_PBS_2   = "ULA_NP_2_counts.txt",
ULA_NTC_PBS_3   = "ULA_NP_3_counts.txt",
ULA_NTC_IFNg_1  = "ULA_NI_1_counts.txt",
ULA_NTC_IFNg_2  = "ULA_NI_2_counts.txt",
ULA_NTC_IFNg_3  = "ULA_NI_3_counts.txt"
)
# -----------------------------
# 3. Read and combine counts
# -----------------------------
count_list <- lapply(count_files, function(f) {
dt <- fread(f, skip = 1)
if (ncol(dt) < 2) stop(paste("File has too few columns:", f))
dt
})
gene_ids <- count_list[[1]]$Geneid
count_list <- lapply(count_list, function(x) x[[ncol(x)]])
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)
# -----------------------------
# 4. Create sample metadata
# -----------------------------
sample_info <- data.frame(
row.names = colnames(count_matrix),
group = factor(c(
rep("Reg_NTC_PBS",3), rep("Reg_NTC_IFNg",3),
rep("ULA_NTC_PBS",3), rep("ULA_NTC_IFNg",3)
))
)
# -----------------------------
# 5. Create DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(
countData = count_matrix,
colData = sample_info,
design = ~ group
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
# -----------------------------
# 6. Extract comparisons
# -----------------------------
# REG comparison
res_reg <- results(dds, contrast = c("group","Reg_NTC_IFNg","Reg_NTC_PBS"))
res_reg <- res_reg[!is.na(res_reg$padj), ]


# ULA comparison
res_ula <- results(dds, contrast = c("group","ULA_NTC_IFNg","ULA_NTC_PBS"))
res_ula <- res_ula[!is.na(res_ula$padj), ]

# Significant DEGs
padj_cut <- 0.05
lfc_cut  <- 1

sig_reg <- rownames(res_reg)[res_reg$padj < padj_cut & abs(res_reg$log2FoldChange) >= lfc_cut]
sig_ula <- rownames(res_ula)[res_ula$padj < padj_cut & abs(res_ula$log2FoldChange) >= lfc_cut]

# Identify unique and shared DEGs
ula_specific <- setdiff(sig_ula, sig_reg)
reg_specific <- setdiff(sig_reg, sig_ula)
shared_genes <- intersect(sig_ula, sig_reg)

cat("DEGs REG-specific:", length(reg_specific), "\n")
cat("DEGs ULA-specific:", length(ula_specific), "\n")
cat("DEGs shared:", length(shared_genes), "\n")

# =========================================
# 6️⃣ Venn diagram
# =========================================
grid.newpage()
venn.plot <- venn.diagram(
  x = list("ULA IFNg vs PBS" = sig_ula, "REG IFNg vs PBS" = sig_reg),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)
grid.draw(venn.plot)
png("DEG_Venn.png", width = 2000, height = 2000, res = 300)
grid.draw(venn.plot)
dev.off()

# =========================================
# 7️⃣ Variance-stabilized counts for heatmap
# =========================================
vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)

# Subset to ULA-specific DEGs
ula_ensembl <- rownames(vst_mat)[rownames(vst_mat) %in% ula_specific]
heatmat <- vst_mat[ula_ensembl, ]

# -----------------------------------------
# Map Ensembl IDs to gene symbols
# -----------------------------------------
gene_map <- bitr(ula_ensembl, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
gene_symbols <- gene_map$SYMBOL
names(gene_symbols) <- gene_map$ENSEMBL
mapped_symbols <- gene_symbols[ula_ensembl]
mapped_symbols[is.na(mapped_symbols)] <- ula_ensembl[is.na(mapped_symbols)]
rownames(heatmat) <- mapped_symbols

# =========================================
# 8️⃣ Column annotations
# =========================================
sample_names <- colnames(heatmat)
treatment <- ifelse(grepl("IFNg", sample_names), "IFNg", "PBS")
group <- ifelse(grepl("ULA", sample_names), "ULA", "REG")
annotation_col <- data.frame(
  Group = group,
  Treatment = treatment
)
rownames(annotation_col) <- sample_names

# =========================================
# 9️⃣ Row annotations: Pathway class
# =========================================
# Define the pathway mapping for each gene (example)
# This must be provided as a named vector: names = gene symbols
# e.g. gene_pathway <- c("Ifit1"="ISG", "Psmb10"="Proteasome / MHC", ...)
# Make sure all ULA-specific genes are included
# Replace this with your actual assignments
gene_pathway <- c(
  Ifi203="ISG", Zbp1="ISG", Isg15="ISG", Gbp9="ISG", Oasl2="ISG",
  Oas3="ISG", Usp18="ISG", Trim30a="ISG", Irf7="ISG", Stat2="ISG",
  Ddx60="ISG", Psmb10="Proteasome / MHC", Irf9="ISG", Irgm1="ISG",
  "9930111J21Rik2"="ISG", Irf1="ISG", Cmpk2="ISG", Apol9b="ISG",
  Ciita="Proteasome / MHC", Rtp4="ISG", Parp14="DNA repair",
  Cd74="Proteasome / MHC", AW112010="ISG", Mpeg1="ISG",
  Cd274="Immune checkpoint", Ifit2="ISG", Ifit3="ISG",
  Ifit1bl1="ISG", Ifit3b="ISG", Ifit1="ISG"
)

# Subset to genes actually in heatmat
gene_pathway <- gene_pathway[rownames(heatmat)]

annotation_row <- data.frame(Pathway = factor(gene_pathway))
rownames(annotation_row) <- names(gene_pathway)

# -----------------------------------------
# Define colors
# -----------------------------------------
pathway_colors <- c(
  "ISG" = "gray",
  "Proteasome / MHC" = "#fee090",
  "DNA repair" = "violet",
  "Immune checkpoint" = "orange"
)
ann_colors <- list(
  Group = c(ULA="purple1", REG="aquamarine"),
  Treatment = c(IFNg="red", PBS="gray"),
  Pathway = pathway_colors
)

# -----------------------------------------
# Order rows by pathway
# -----------------------------------------
ordered_genes <- rownames(annotation_row)[order(annotation_row$Pathway)]
heatmat_ordered <- heatmat[ordered_genes, ]
annotation_row <- annotation_row[ordered_genes, , drop = FALSE]

# =========================================
# 10️⃣ Draw annotated heatmap
# =========================================
pheatmap(
  heatmat_ordered,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  fontsize_row = 8,
  main = "ULA-specific IFNg DEGs grouped by functional pathways"
)

# Export the final VST matrix used in the heatmap (with gene symbol rownames)
write.csv(heatmat_ordered, "Fig_4D_Heatmap_VST_Matrix.csv", row.names = TRUE)
