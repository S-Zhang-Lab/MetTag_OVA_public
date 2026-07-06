# ============================
# Curated receptor analysis: Reg PBS vs ULA IFNg
# ============================
# === Curated cytokine receptor heatmap: Reg PBS vs ULA IFNg ===
library(DESeq2)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)
library(dplyr)
library(pheatmap)

# -----------------------------
# 1️⃣ Define count files
# -----------------------------
count_files <- list(
  Reg_NTC_PBS_1 = "Reg_NP_1_counts.txt",   
  Reg_NTC_PBS_2 = "Reg_NP_2_counts.txt",
  Reg_NTC_PBS_3 = "Reg_NP_3_counts.txt",
  ULA_NTC_IFNg_1  = "ULA_NI_1_counts.txt",  
  ULA_NTC_IFNg_2  = "ULA_NI_2_counts.txt",
  ULA_NTC_IFNg_3  = "ULA_NI_3_counts.txt"
)

# -----------------------------
# 2️⃣ Read counts
# -----------------------------
count_list <- lapply(count_files, function(f) fread(f, skip = 1)[[ncol(fread(f, skip = 1))]])
gene_ids <- fread(count_files[[1]], skip=1)$Geneid
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)

# -----------------------------
# 3️⃣ Sample metadata
# -----------------------------
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c(rep("Reg_PBS",3), rep("ULA_IFNg",3)))
)
sample_info$condition <- relevel(sample_info$condition, ref = "Reg_PBS")

# -----------------------------
# 4️⃣ DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # pre-filter low counts
dds <- DESeq(dds)

# -----------------------------
# 5️⃣ DESeq2 results
# -----------------------------
resLFC <- results(dds)  # this creates the object you were missing
resLFC_df <- as.data.frame(resLFC)
resLFC_df$ENSEMBL <- rownames(resLFC_df)

# -----------------------------
# 6️⃣ Map Ensembl -> SYMBOL
# -----------------------------
ensembl2symbol <- bitr(resLFC_df$ENSEMBL,
                       fromType="ENSEMBL",
                       toType="SYMBOL",
                       OrgDb=org.Mm.eg.db)
vst_mat <- assay(vst(dds, blind = FALSE))
vst_sym <- vst_mat
rownames(vst_sym) <- ensembl2symbol$SYMBOL[match(rownames(vst_mat), ensembl2symbol$ENSEMBL)]

# -----------------------------
# 7️⃣ Curated receptor list
# -----------------------------
custom_receptors <- c(
  "Tgfbr1","Tgfbr2","Tgfbr3",
  "Tnfrsf1a","Tnfrsf1b",
  "Il1r1","Il1r2","Il1rap",
  "Il10ra","Il10rb"
)

# -----------------------------
# 8️⃣ Subset VST to receptors present
# -----------------------------
receptors_present <- intersect(custom_receptors, rownames(vst_sym))
if(length(receptors_present) == 0) stop("No curated receptors found in VST matrix!")
cat("Receptors present:\n")
print(receptors_present)

vst_receptors <- vst_sym[receptors_present, , drop = FALSE]

# -----------------------------
# 9️⃣ Row-scale for heatmap
# -----------------------------
vst_receptors_scaled <- t(scale(t(vst_receptors)))

# -----------------------------
# 🔟 Annotation & Heatmap
# -----------------------------
annotation_col <- data.frame(Condition = colData(dds)$condition)
rownames(annotation_col) <- colnames(vst_receptors_scaled)

pheatmap(vst_receptors_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         main = "Cytokine receptor expression: Reg PBS vs ULA IFNg",
         fontsize_row = 10,
         show_rownames = TRUE)

# --- 1️⃣ Extract DESeq2 results for your comparison of interest ---
# ULA_IFNg vs Reg_PBS
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)

# --- DESeq2 results for your contrast ---
res <- results(dds, contrast = c("condition", "ULA_IFNg", "Reg_PBS"))

# --- Convert to data frame ---
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column(var = "ENSEMBL")

# --- Map Ensembl IDs to SYMBOLs ---
symbol_map <- bitr(res_df$ENSEMBL, fromType = "ENSEMBL",
                   toType = "SYMBOL", OrgDb = org.Mm.eg.db)

res_df <- left_join(res_df, symbol_map, by = "ENSEMBL")

# --- Curated receptor list ---
custom_receptors <- c("Tgfbr1","Tgfbr2","Tgfbr3",
                      "Tnfrsf1a","Tnfrsf1b",
                      "Il1r1","Il1r2","Il1rap",
                      "Il10ra","Il10rb")

# --- Subset DESeq2 results to receptors ---
res_receptors_df <- res_df %>%
  filter(SYMBOL %in% custom_receptors)

# --- Extract VST ---
vst_mat <- assay(vst(dds, blind = FALSE))
# Map VST rownames if necessary (Ensembl → SYMBOL)
vst_df <- as.data.frame(vst_mat) %>%
  tibble::rownames_to_column(var = "ENSEMBL")
vst_df <- left_join(vst_df, symbol_map, by = "ENSEMBL")

# --- Keep only receptors ---
library(dplyr)
library(tibble)

vst_receptors <- vst_df %>%
  filter(SYMBOL %in% custom_receptors) %>%
  as_tibble() %>%           # convert to tibble
  dplyr::select(SYMBOL, dplyr::everything(), -ENSEMBL) %>%  # explicit dplyr::select
  column_to_rownames(var = "SYMBOL")


# --- Convert to long format for summary stats ---
expr_long <- vst_receptors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  pivot_longer(cols = -SYMBOL, names_to = "Sample", values_to = "VST") %>%
  left_join(as.data.frame(colData(dds)) %>% tibble::rownames_to_column("Sample"), by = "Sample")

# --- Compute per-gene summary stats ---
gene_stats <- expr_long %>%
  group_by(SYMBOL) %>%
  summarize(
    mean_Reg_PBS = mean(VST[condition == "Reg_PBS"], na.rm = TRUE),
    mean_ULA_IFNg = mean(VST[condition == "ULA_IFNg"], na.rm = TRUE),
    med_Reg_PBS = median(VST[condition == "Reg_PBS"], na.rm = TRUE),
    med_ULA_IFNg = median(VST[condition == "ULA_IFNg"], na.rm = TRUE),
    wilcox_p = tryCatch(wilcox.test(VST ~ condition, data = cur_data_all())$p.value,
                        error = function(e) NA_real_),
    n_Reg_PBS = sum(condition == "Reg_PBS"),
    n_ULA_IFNg = sum(condition == "ULA_IFNg"),
    .groups = "drop"
  ) %>%
  mutate(wilcox_p_adj = p.adjust(wilcox_p, method = "BH"))

# --- Combine DESeq2 stats with summary stats ---
library(dplyr)
library(tibble)

# --- Combine DESeq2 stats with summary stats ---
final_table <- dplyr::left_join(
  gene_stats,
  res_receptors_df %>%
    dplyr::as_tibble() %>%
    dplyr::select(SYMBOL, log2FoldChange, lfcSE, stat, pvalue, padj),
  by = "SYMBOL"
)

# --- Export ---
write.csv(final_table, "curated_receptor_ULAIFNg_vs_RegPBS.csv", row.names = FALSE)

# ============================
# Curated receptor analysis: 4 NTC groups (Reg vs ULA, PBS vs IFNg)
# ============================

# --- Load libraries ---
library(DESeq2)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)
library(dplyr)
library(pheatmap)
library(tidyr)

setwd("/Users/s438978/Desktop/RData/RNA-seq")

# -----------------------------
# 1️⃣ Define count files (4 NTC groups only)
# -----------------------------
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
# 2️⃣ Read count matrices
# -----------------------------
count_list <- lapply(count_files, function(f) fread(f, skip = 1)[[ncol(fread(f, skip = 1))]])
gene_ids <- fread(count_files[[1]], skip = 1)$Geneid
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- names(count_files)

# -----------------------------
# 3️⃣ Build sample metadata
# -----------------------------
sample_info <- data.frame(
  Sample = colnames(count_matrix)
) %>%
  mutate(
    Line = ifelse(grepl("Reg", Sample), "Reg", "ULA"),
    Genotype = "NTC",
    Treatment = ifelse(grepl("IFNg", Sample), "IFNg", "PBS"),
    condition = paste(Line, Genotype, Treatment, sep = "_")
  )

rownames(sample_info) <- sample_info$Sample

# -----------------------------
# 4️⃣ DESeq2 setup
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# -----------------------------
# 5️⃣ VST transformation & mapping
# -----------------------------
vst_mat <- assay(vst(dds, blind = FALSE))
ensembl_ids <- rownames(vst_mat)
symbol_map <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
vst_sym <- vst_mat
rownames(vst_sym) <- symbol_map$SYMBOL[match(rownames(vst_mat), symbol_map$ENSEMBL)]
vst_sym <- vst_sym[!is.na(rownames(vst_sym)), ]

# -----------------------------
# 6️⃣ Updated curated receptor list
# -----------------------------
custom_receptors <- c(
  "Tgfbr1","Tgfbr2","Tgfbr3","Pdgfra",
  "Tnfrsf21","Tnfrsf1b","Ltbr","Tnfrsf1a",
  "Itgb8","Eng","Acvrl1","App","Adrb2",
  "Il12rb2","Il1r1","Il1r2","Traf2",
  "Il10ra","Il10rb"
)

# -----------------------------
# 7️⃣ Subset to receptors & row-scale
# -----------------------------
receptors_present <- intersect(custom_receptors, rownames(vst_sym))
cat("Receptors present in dataset:\n")
print(receptors_present)

vst_receptors <- vst_sym[receptors_present, , drop = FALSE]
vst_receptors_scaled <- t(scale(t(vst_receptors)))

# -----------------------------
# 8️⃣ Order samples by Treatment then Line
# -----------------------------
sample_order <- sample_info %>%
  arrange(Treatment, Line) %>%
  pull(Sample)

vst_receptors_scaled <- vst_receptors_scaled[, sample_order]
annotation_col <- sample_info[sample_order, c("Line","Treatment")]

# -----------------------------
# 9️⃣ Heatmap: grouped by treatment and condition
# -----------------------------
pheatmap(vst_receptors_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,     # no clustering, keep order
         annotation_col = annotation_col,
         main = "Curated receptor expression (NTC groups only)",
         fontsize_row = 10,
         show_rownames = TRUE,
         angle_col = 45)

# -----------------------------
# 🔟 Export receptor expression table
# -----------------------------
write.csv(as.data.frame(vst_receptors),
          "curated_receptor_expression_NTC_groups.csv",
          row.names = TRUE)

### check St6gal1
library(org.Mm.eg.db)

st6_info <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = "St6gal1",
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)

st6_id <- st6_info$ENSEMBL[1]
st6_id

norm_counts <- counts(dds, normalized = TRUE)

# Extract expression for St6gal1
st6_expr <- norm_counts[st6_id, ]

df_st6 <- data.frame(
  expression = as.numeric(st6_expr),
  condition = sample_info[colnames(norm_counts), "condition"],
  Line = sample_info[colnames(norm_counts), "Line"],
  Treatment = sample_info[colnames(norm_counts), "Treatment"]
)

library(ggplot2)

df_st6$condition <- factor(
  df_st6$condition,
  levels = c(
    "Reg_NTC_PBS",
    "Reg_NTC_IFNg",
    "ULA_NTC_PBS",
    "ULA_NTC_IFNg"
  )
)


ggplot(df_st6, aes(x = condition, y = expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.15), size = 3) +
  scale_y_log10() +
  labs(title = "St6gal1 Expression Across Four NTC Conditions",
       x = "Condition",
       y = "Normalized Expression (log10)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

group1 <- df_st6 %>% filter(condition == "Reg_NTC_PBS") %>% pull(expression)
group2 <- df_st6 %>% filter(condition == "ULA_NTC_PBS") %>% pull(expression)

wilcox_res <- wilcox.test(group1, group2)
wilcox_res

t.test(group1, group2)
# Welch Two Sample t-test
#data:  group1 and group2
#t = -13.353, df = 2.2841, p-value = 0.003238
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -6336.026 -3512.451
#sample estimates:
#  mean of x mean of y 
#3010.288  7934.526 

group1 <- df_st6 %>% filter(condition == "ULA_NTC_PBS") %>% pull(expression)
group2 <- df_st6 %>% filter(condition == "ULA_NTC_IFNg") %>% pull(expression)

# NS


