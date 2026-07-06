################################################################################
# Script: 12_IFNg_Scoring.R
# Purpose: Evaluation of IFN responses, clonal time point divergence (TP-1 vs TP-6),
#          and downstream stress/survival programs in Ovarian Cancer Metastases.
################################################################################

# ==============================================================================
# PART 1: INITIALIZATION & SIGNATURE SETUP
# ==============================================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(pheatmap)
library(clusterProfiler)
library(writexl)

set.seed(42)
setwd("C:/Users/Zhang Lab/Desktop/Emma/R01_and_2_w_LARRY")

# 1. Load the Seurat Object
Mets <- readRDS("Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")

# 2. Define IFN Signatures
# Extract genes for the INTERFERON_ALPHA_RESPONSE set
hallmark_gmt <- read.gmt("mh.all.v2024.1.Mm.symbols.gmt.txt")
ifna_genes <- hallmark_gmt %>%
  filter(term == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene)

# Define Canonical IFN-gamma genes
canonical_ifng_genes <- c(
  "Irf1", "Irf8", "Stat1", "Stat2", "Socs1", "Ciita", "Cxcl9", "Cxcl10", 
  "Cxcl11", "Cxcl13", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Psme1", 
  "Psme2", "H2-K1", "H2-D1", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "B2m",
  "Gbp2b", "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Gbp6", "Gbp7", "Isg15", "Ifi30", 
  "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oas2", "Oas3", "Oasl1", 
  "Rsad2", "Irf7", "Lgals3bp", "Trim21", "Nos2", "Casp1", "Casp8", "Tnf", 
  "Cd274", "Ifngr1", "Ifngr2"
)

# 3. Add Module Scores to Metadata
Mets <- AddModuleScore(Mets, features = list(canonical_ifng_genes), name = "IFNG_score")
Mets <- AddModuleScore(Mets, features = list(ifna_genes), name = "IFNa_score")


# ==============================================================================
# PART 2: GLOBAL COMPARISON (AscMet vs OmMet)
# ==============================================================================
# 1. Downsample to match sample sizes
df_meta <- Mets@meta.data
ascmet_size <- nrow(df_meta[df_meta$Combined.HTO_group == "asc:1_AscMet", ])

df_global_downsampled <- df_meta %>%
  group_by(Combined.HTO_group) %>%
  group_modify(~ {
    if (.y$Combined.HTO_group == "asc:1_AscMet") return(.x)
    else sample_n(.x, size = ascmet_size)
  }) %>%
  ungroup() %>%
  mutate(Combined.HTO_group = factor(Combined.HTO_group))

# 2. Calculate Statistics
pval_IFNG_global <- wilcox.test(IFNG_score1 ~ Combined.HTO_group, data = df_global_downsampled)$p.value
pval_IFNa_global <- wilcox.test(IFNa_score1 ~ Combined.HTO_group, data = df_global_downsampled)$p.value

# 3. Plotting Function
plot_global_score <- function(data, score_col, y_label, pval) {
  ggplot(data, aes(x = Combined.HTO_group, y = .data[[score_col]], fill = Combined.HTO_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    annotate("text", x = 1.5, y = max(data[[score_col]], na.rm = TRUE) * 1.05, 
             label = paste0("p = ", signif(pval, 3)), size = 5) +
    labs(y = y_label, x = "Sample Group") +
    theme_minimal() + theme(legend.position = "none")
}

p1 <- plot_global_score(df_global_downsampled, "IFNG_score1", "IFNG Signature Score", pval_IFNG_global)
p2 <- plot_global_score(df_global_downsampled, "IFNa_score1", "IFNa Signature Score", pval_IFNa_global)
p1 + p2

# 4. Export Source Data
global_source_data <- list(
  "FigX_Raw_Scores" = df_global_downsampled %>% select(Combined.HTO_group, IFNG_score1, IFNa_score1) %>% rownames_to_column("Cell_Barcode"),
  "FigX_Statistics" = data.frame(
    Signature = c("IFNG_score1", "IFNa_score1"),
    N_Cells_Per_Group = c(ascmet_size, ascmet_size),
    P_Value = c(pval_IFNG_global, pval_IFNa_global)
  )
)
write_xlsx(global_source_data, "Source_Data_Global_IFN_Signatures.xlsx")



# PART 3: TRAJECTORY TIME POINT ANALYSIS (TP-1 through TP-6)
# ==============================================================================
# 1. Map All 6 Time Points based on Omental Barcodes
df_tp_all <- Mets@meta.data %>%
  filter(om_LibID_first_bc %in% c("R1_BC_ID6", "R1_BC_ID5", "R1_BC_ID4", "R1_BC_ID3", "R1_BC_ID2", "R1_BC_ID1")) %>%
  mutate(
    Time_Point = case_when(
      om_LibID_first_bc == "R1_BC_ID6" ~ "TP-1 (BC ID6)",
      om_LibID_first_bc == "R1_BC_ID5" ~ "TP-2 (BC ID5)",
      om_LibID_first_bc == "R1_BC_ID4" ~ "TP-3 (BC ID4)",
      om_LibID_first_bc == "R1_BC_ID3" ~ "TP-4 (BC ID3)",
      om_LibID_first_bc == "R1_BC_ID2" ~ "TP-5 (BC ID2)",
      om_LibID_first_bc == "R1_BC_ID1" ~ "TP-6 (BC ID1)"
    )
  ) %>%
  # Lock the factor order chronologically
  mutate(Time_Point = factor(Time_Point, levels = c(
    "TP-1 (BC ID6)", "TP-2 (BC ID5)", "TP-3 (BC ID4)", 
    "TP-4 (BC ID3)", "TP-5 (BC ID2)", "TP-6 (BC ID1)"
  )))

# 2. Calculate Statistics (Strictly TP-1 vs TP-6)
# Filter just the terminal ends for the math
df_tp_ends <- df_tp_all %>% filter(Time_Point %in% c("TP-1 (BC ID6)", "TP-6 (BC ID1)"))

p_ifna_tp <- wilcox.test(IFNa_score1 ~ Time_Point, data = df_tp_ends)$p.value
p_ifng_tp <- wilcox.test(IFNG_score1 ~ Time_Point, data = df_tp_ends)$p.value

cat("Wilcoxon p-value (IFNa: TP-1 vs TP-6):", signif(p_ifna_tp, 3), "\n")
cat("Wilcoxon p-value (IFNg: TP-1 vs TP-6):", signif(p_ifng_tp, 3), "\n")

# 3. Export Source Data (All 6 TPs for the plot, 2 TPs for the stats)
tp_source_data <- list(
  "FigX_TimePoint_Violins" = df_tp_all %>% 
    rownames_to_column("Cell_Barcode") %>% 
    select(Cell_Barcode, Time_Point, IFNa_score1, IFNG_score1) %>% 
    arrange(Time_Point),
  
  "FigX_TimePoint_Stats" = data.frame(
    Feature = c("IFNa_score1", "IFNG_score1"),
    Comparison = c("TP-1 vs TP-6", "TP-1 vs TP-6"),
    P_Value = c(p_ifna_tp, p_ifng_tp)
  )
)
write_xlsx(tp_source_data, "Source_Data_TP1_through_TP6_FINAL.xlsx")
# ==============================================================================
# PART 4: IFNg RESPONDER GATING & DOWNSTREAM PROGRAMS
# ==============================================================================
# 1. Manual Gating via FeatureScatter (Interactive)
Idents(Mets) <- "seurat_clusters"
plotActbxIrf1 <- FeatureScatter(Mets, feature1 = "Actb", feature2 = "Irf1")
cat("Please select IFNg-low cells, then IFNg-high cells in the interactive prompt...\n")
Mets <- CellSelector(plot = plotActbxIrf1, object = Mets, ident = "IFNg-low")
IFNg_low_cells <- WhichCells(Mets, idents = "IFNg-low")

Mets <- CellSelector(plot = plotActbxIrf1, object = Mets, ident = "IFNg-high")
IFNg_high_cells <- WhichCells(Mets, idents = "IFNg-high")

# Update object metadata with manual labels
Mets$IFNg_status <- "Mid"
Mets$IFNg_status[IFNg_low_cells] <- "IFNg-low"
Mets$IFNg_status[IFNg_high_cells] <- "IFNg-high"

# 2. Define Downstream Programs
autophagy_genes   <- c("Atg3", "Atg5", "Atg7", "Atg12", "Atg16l1", "Becn1", "Becn2", "Map1lc3b", "Sqstm1", "Ulk1", "Ulk2", "Ambra1", "Wipi1", "Gabarapl1", "Gabarapl2", "Gabarap", "Rragd", "Tsc1", "Tsc2", "Mtor", "Bnip3", "Bnip3l", "Lamp1", "Lamp2", "Pik3c3", "Rraga", "Rragb", "Tfeb", "Foxo3")
prosurvival_genes <- c("Bcl2", "Bcl2l1", "Mcl1", "Xiap", "Birc2", "Birc3", "Cflar", "Tnfrsf10b", "Birc5", "Bcl2a1", "Stat3", "Akt1")
dormancy_genes    <- c("Cdkn1a", "Cdkn2a", "Gadd45a", "Gadd45g", "Ddit3", "Fos", "Junb")
stress_genes      <- c("Hspa1a", "Hspa1b", "Atf4", "Atf3", "Hspb1", "Trib3", "Xbp1", "Sod1", "Sod2", "Gpx1", "Nfe2l2", "Hif1a")

# 3. Add Scores
Mets <- AddModuleScore(Mets, features = list(autophagy_genes), name = "AutophagyScore")
Mets <- AddModuleScore(Mets, features = list(prosurvival_genes), name = "ProSurvivalScore")
Mets <- AddModuleScore(Mets, features = list(dormancy_genes), name = "DormancyScore")
Mets <- AddModuleScore(Mets, features = list(stress_genes), name = "StressScore")

# 4. Extract and Compare High vs Low Responders
responder_cells <- c(IFNg_low_cells, IFNg_high_cells)
df_programs <- FetchData(Mets, vars = c("AutophagyScore1", "ProSurvivalScore1", "DormancyScore1", "StressScore1", "IFNg_status"), cells = responder_cells)

p_auto <- wilcox.test(AutophagyScore1 ~ IFNg_status, data = df_programs)$p.value
p_surv <- wilcox.test(ProSurvivalScore1 ~ IFNg_status, data = df_programs)$p.value
p_dorm <- wilcox.test(DormancyScore1 ~ IFNg_status, data = df_programs)$p.value
p_stre <- wilcox.test(StressScore1 ~ IFNg_status, data = df_programs)$p.value

# Export Source Data
program_source_data <- list(
  "FigX_Program_Scores" = df_programs %>% rownames_to_column("Cell_Barcode") %>% arrange(IFNg_status),
  "FigX_Program_Stats" = data.frame(
    Program = c("Autophagy", "Apoptosis/Survival", "Dormancy", "Stress"),
    P_Value = c(p_auto, p_surv, p_dorm, p_stre)
  )
)
write_xlsx(program_source_data, "Source_Data_IFNg_Downstream_Programs.xlsx")


# ==============================================================================
# PART 5: DIFFERENTIAL EXPRESSION & TARGET HEATMAPS
# ==============================================================================
# 1. Target Gene Expression (Irf1 High vs Low)
Idents(Mets) <- "IFNg_status"
subset_Mets <- subset(Mets, idents = c("IFNg-high", "IFNg-low"))
DefaultAssay(subset_Mets) <- "RNA"
subset_Mets <- NormalizeData(subset_Mets)

genes_to_test <- c("Psmb10", "Psmb9", "Parp12", "Parp14", "Parp9", "Cd274", "Mcl1", "Bcl2", "Socs3")
de_results_targets <- FindMarkers(subset_Mets, ident.1 = "IFNg-high", ident.2 = "IFNg-low", features = genes_to_test, test.use = "wilcox")

# 2. Sample-Level Heatmaps (Z-Scores)
# Filter singlets
Mets_singlet <- subset(Mets, subset = om_HTO_classification.global == "Singlet" | Asc_HTO_classification.global == "Singlet")
Mets_singlet$HTO_ID <- ifelse(!is.na(Mets_singlet$om_hash.ID), Mets_singlet$om_hash.ID, Mets_singlet$Asc_hash.ID)
hto_ids <- factor(Mets_singlet$HTO_ID)

expr_sct <- GetAssayData(Mets_singlet, assay = "SCT", slot = "data")
valid_ifna <- intersect(ifna_genes, rownames(expr_sct))
valid_ifng <- intersect(canonical_ifng_genes, rownames(expr_sct))

average_expr_by_sample <- function(genes, expr_mat, sample_ids) {
  sapply(levels(sample_ids), function(s) rowMeans(expr_mat[genes, sample_ids == s, drop = FALSE]))
}

# Generate scaled matrices
scaled_ifna <- t(scale(t(average_expr_by_sample(valid_ifna, expr_sct, hto_ids))))
scaled_ifng <- t(scale(t(average_expr_by_sample(valid_ifng, expr_sct, hto_ids))))

# Plot Heatmaps
pheatmap(scaled_ifna, cluster_rows = TRUE, cluster_cols = TRUE, main = "IFN-alpha Genes", fontsize_row = 8)
pheatmap(scaled_ifng, cluster_rows = TRUE, cluster_cols = TRUE, main = "IFN-gamma Genes", fontsize_row = 8)

# 3. Final Export for Targets and Heatmaps
target_source_data <- list(
  "FigX_Irf1_DEGs" = de_results_targets %>% rownames_to_column("Gene"),
  "FigX_IFNa_Heatmap" = as.data.frame(scaled_ifna) %>% rownames_to_column("Gene"),
  "FigX_IFNg_Heatmap" = as.data.frame(scaled_ifng) %>% rownames_to_column("Gene")
)
write_xlsx(target_source_data, "Source_Data_Target_Genes_and_Heatmaps.xlsx")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
