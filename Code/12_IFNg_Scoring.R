### EA additional analysis of SZ's Mets object ###

# load required packages
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(scales)
library(plotly)
library(pheatmap)
library(viridis) 
library(clusterProfiler)

# set working directory
setwd("/Users/s438978/Desktop/RData/MetTag_U6_LARRY_BC_OVA/scVelo2_SZ")

# load rds saved from Seurat
Mets <- readRDS("Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")
DimPlot(Mets)

# === Examine IFNa and IFNg genesets ================
# Read hallmark gene set
hallmark_gmt <- read.gmt("/Users/s438978/Desktop/RData/MetTag_U6_LARRY_BC_OVA/scVelo2_SZ/mh.all.v2024.1.Mm.symbols.gmt.txt")

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

Mets<- AddModuleScore(Mets, features = list(canonical_ifng_genes), name = "IFNG_score")
VlnPlot(Mets, features = "IFNG_score1", group.by = "Combined.HTO_group")

Mets<- AddModuleScore(Mets, features = list(ifna_genes), name = "IFNa_score")
VlnPlot(Mets, features = "IFNa_score1", group.by = "Combined.HTO_group")

library(dplyr)
library(ggplot2)
library(ggpubr)  # for easy p-value annotation

# Example dataframe: df with columns Combined.HTO_group, Score1, Score2
# Replace `df` with your actual dataframe name

# Step 1: Downsample larger groups to match the size of asc:1_AscMet
df <- Mets@meta.data
ascmet_cells <- df %>% filter(Combined.HTO_group == "asc:1_AscMet")
ascmet_size <- nrow(ascmet_cells)

# Downsample other groups to ascmet_size each
df_downsampled <- df %>%
  group_by(Combined.HTO_group) %>%
  group_modify(~ {
    if (.y$Combined.HTO_group == "asc:1_AscMet") {
      return(.x)
    } else {
      sample_n(.x, size = ascmet_size)
    }
  }) %>%
  ungroup()

# Step 2: Plot the two scores with p-values between groups

# Convert Combined.HTO_group to factor if needed
df_downsampled$Combined.HTO_group <- factor(df_downsampled$Combined.HTO_group)

# Define a function to plot a score with p-value between groups
plot_score_with_pval <- function(data, score_col, y_label) {
  ggplot(data, aes(x = Combined.HTO_group, y = .data[[score_col]], fill = Combined.HTO_group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",       # shows formatted p-value, e.g. 1.23e-10
      label.x.npc = "center",   # position label at center on x-axis
      label.y.npc = 0.95,       # position label near top (95% of y-axis range)
      hide.ns = TRUE            # hides p-value if not significant
    ) +
    labs(y = y_label, x = "Sample Group") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Plot Score1 and Score2
p1 <- plot_score_with_pval(df_downsampled, "IFNG_score1", "IFNG_score1")
p2 <- plot_score_with_pval(df_downsampled, "IFNa_score1", "IFNa_score1")


# Display plots side by side
library(patchwork)
p1 + p2

### since p values are really small, add them manually
summary(df_downsampled$IFNG_score1)
summary(df_downsampled$IFNa_score1)
pval_IFNG <- wilcox.test(IFNG_score1 ~ Combined.HTO_group, data = df_downsampled)$p.value
pval_IFNa <- wilcox.test(IFNa_score1 ~ Combined.HTO_group, data = df_downsampled)$p.value

plot_score_with_manual_pval <- function(data, score_col, y_label, pval) {
  ggplot(data, aes(x = Combined.HTO_group, y = .data[[score_col]], fill = Combined.HTO_group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    annotate(
      "text", x = 1.5, y = max(data[[score_col]], na.rm = TRUE) * 1.05, 
      label = paste0("p = ", format(pval, scientific = TRUE, digits = 3)),
      size = 5
    ) +
    labs(y = y_label, x = "Sample Group") +
    theme_minimal() +
    theme(legend.position = "none")
}

p1 <- plot_score_with_manual_pval(df_downsampled, "IFNG_score1", "IFNG Signature Score", pval_IFNG)
p2 <- plot_score_with_manual_pval(df_downsampled, "IFNa_score1", "IFNa Signature Score", pval_IFNa)

library(patchwork)
p1 + p2

# ==== compare DEGs between cells of Mets@meta.data$om_LibID_first_bc == R1_BC_ID6 and Mets@meta.data$Asc_LibID_first_bc == R1_BC_ID6  cells. 
# Create a new identity class in the Seurat object to define these two populations
# Initialize with NA
Mets@meta.data$LibID_bc_group <- NA

# Assign Om group: R1_BC_ID6 or R1_BC_ID5 in om_LibID_first_bc
Mets@meta.data$LibID_bc_group[Mets@meta.data$om_LibID_first_bc %in% c("R1_BC_ID6", "R1_BC_ID5")] <- "Om_R1_BC_ID6_or_5"

# Assign Asc group: R1_BC_ID6 or R1_BC_ID5 in Asc_LibID_first_bc
Mets@meta.data$LibID_bc_group[Mets@meta.data$Asc_LibID_first_bc %in% c("R1_BC_ID6", "R1_BC_ID5")] <- "Asc_R1_BC_ID6_or_5"

# Check counts
table(Mets$LibID_bc_group)

# Set identities based on this grouping
Idents(Mets) <- "LibID_bc_group"

# Subset to keep only cells with non-NA groups
cells_to_keep <- rownames(Mets@meta.data)[!is.na(Mets@meta.data$LibID_bc_group)]
Mets_filtered <- subset(Mets, cells = cells_to_keep)

# Set Idents again for filtered object
Idents(Mets_filtered) <- "LibID_bc_group"

#--------------------
# 1. Add the grouping info to metadata (assuming group_vector is ready)
group_df <- data.frame(LibID_bc_group = group_vector)
rownames(group_df) <- names(group_vector)

# Add/update metadata in the Seurat object
Mets@meta.data <- cbind(Mets@meta.data, group_df)

# 2. Define groups to keep (TP-1 through TP-6)
groups_to_keep <- paste0("TP-", 1:6)

# 3. Subset to keep only these cells
cells_to_keep <- WhichCells(Mets, expression = LibID_bc_group %in% groups_to_keep)
Mets_filtered <- subset(Mets, cells = cells_to_keep)

# 4. Create a clean factor for plotting
Mets_filtered$TP_label <- factor(Mets_filtered$LibID_bc_group, levels = groups_to_keep)

# 5. Set Idents to TP_label for plotting convenience
Idents(Mets_filtered) <- Mets_filtered$TP_label

# 6. Plot IFNa and IFNg signature scores side by side with rotated x labels
VlnPlot(Mets_filtered,
        features = c("IFNa_score1", "IFNG_score1"),
        pt.size = 0.2, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 7. Run Wilcoxon test comparing TP-1 vs TP-6 for IFNa_score1 and IFNG_score1
library(dplyr)
df <- Mets_filtered@meta.data %>%
  filter(TP_label %in% c("TP-1", "TP-6"))

p_ifna <- wilcox.test(IFNa_score1 ~ TP_label, data = df)$p.value
cat("Wilcoxon p-value (IFNa_score1: TP-1 vs TP-6):", format(p_ifna, scientific = TRUE), "\n")
#Wilcoxon p-value (IFNa_score1: TP-1 vs TP-6): 1.276725e-01

p_ifng <- wilcox.test(IFNG_score1 ~ TP_label, data = df)$p.value
cat("Wilcoxon p-value (IFNG_score1: TP-1 vs TP-6):", format(p_ifng, scientific = TRUE), "\n")
# Wilcoxon p-value (IFNG_score1: TP-1 vs TP-6): 4.502236e-03

# Subset metadata
df <- Mets_filtered@meta.data %>% filter(TP_label %in% c("TP-1", "TP-6"))

# Define features
features_to_test <- c("IFNa_score1", "IFNG_score1")

# Initialize results list
results_list <- lapply(features_to_test, function(feature) {
  test_result <- wilcox.test(as.formula(paste(feature, "~ TP_label")), data = df)
  
  # Calculate log fold change (mean TP-1 - mean TP-6)
  tp1_mean <- mean(df[[feature]][df$TP_label == "TP-1"], na.rm = TRUE)
  tp6_mean <- mean(df[[feature]][df$TP_label == "TP-6"], na.rm = TRUE)
  log_fc <- log2(tp1_mean + 1) - log2(tp6_mean + 1)  # log fold change with +1 to avoid log(0)
  
  data.frame(
    Feature = feature,
    P_Value = test_result$p.value,
    Log2_FC_TP1_vs_TP6 = log_fc
  )
})

# Combine into one data frame
results_df <- do.call(rbind, results_list)

# Export to CSV
write.csv(results_df, "wilcoxon_results_TP1_vs_TP6.csv", row.names = FALSE)

Idents(Mets) <- "seurat_clusters"
VlnPlot(Mets, features = c("Tgfbr1", "Tgfbr2", "Tgfbr3", "Pdgfra", "Itgb8", "Eng", "Acvrl1", "App", "Itgav", "Itgb1", "Itgb3", "Sdc2"))

### 07_GATING out IFNg responsive cells in Mets ###

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
Idents(Mets) <- "seurat_clusters"

# 2. Create your scatter plot for gating on Stat1 vs Irf1
plotActbxIrf1 <- FeatureScatter(Mets, feature1 = "Actb", feature2 = "Irf1")

# 3. Use CellSelector to manually define IFNg-low first
Mets <- CellSelector(plot = plotActbxIrf1, object = Mets, ident = "IFNg-low")

# 4. Save those cells before moving on!
IFNg_low_cells <- WhichCells(Mets, idents = "IFNg-low")

# 5. Now select the IFNg-high cells
Mets <- CellSelector(plot = plotActbxIrf1, object = Mets, ident = "IFNg-high")

# 6. Save those cells too
IFNg_high_cells <- WhichCells(Mets, idents = "IFNg-high")

# 7. Now you have both cell groups!
autophagy_genes <- c("Atg3", "Atg5", "Atg7", "Atg12", "Atg16l1", "Becn1", "Becn2", "Map1lc3b", 
                     "Sqstm1", "Ulk1", "Ulk2", "Ambra1", "Wipi1", "Gabarapl1", "Gabarapl2", 
                     "Gabarap", "Rragd", "Tsc1", "Tsc2", "Mtor", "Bnip3", "Bnip3l", 
                     "Lamp1", "Lamp2", "Pik3c3", "Rraga", "Rragb", "Tfeb", "Foxo3")

Mets <- AddModuleScore(Mets, features = list(autophagy_genes), name = "AutophagyScore")

# Extract autophagy scores
IFNlow_score <- FetchData(Mets, vars = "AutophagyScore1", cells = IFNg_low_cells)
IFNhigh_score <- FetchData(Mets, vars = "AutophagyScore1", cells = IFNg_high_cells)

# Combine into a data frame for plotting
autophagy_df <- data.frame(
  Score = c(IFNlow_score$AutophagyScore1, IFNhigh_score$AutophagyScore1),
  Group = c(rep("IFNg-low", length(IFNlow_score$AutophagyScore1)),
            rep("IFNg-high", length(IFNhigh_score$AutophagyScore1)))
)

# Make a simple violin plot
library(ggplot2)
ggplot(autophagy_df, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white") +
  theme_minimal() +
  ylab("Autophagy Score") +
  xlab("") +
  ggtitle("Autophagy in IFNg Responding vs Non-Responding Cancer Cells")

# ---------------------------
# 1. Define gene lists
# ---------------------------

# Apoptosis Resistance Program
prosurvival_genes <- c("Bcl2", "Bcl2l1", "Mcl1", "Xiap", "Birc2", "Birc3", "Cflar", "Tnfrsf10b", "Birc5", "Bcl2a1", "Stat3", "Akt1")

# Dormancy/Quiescence Program
dormancy_genes <- c("Cdkn1a", "Cdkn2a", "Gadd45a", "Gadd45g", "Ddit3", "Fos", "Junb")

# Stress Tolerance Program
stress_genes <- c("Hspa1a", "Hspa1b", "Atf4", "Atf3", "Hspb1", "Trib3", "Xbp1", "Sod1", "Sod2", "Gpx1", "Nfe2l2", "Hif1a")

# ---------------------------
# 2. Add module scores
# ---------------------------

Mets <- AddModuleScore(Mets, features = list(prosurvival_genes), name = "ProSurvivalScore")
Mets <- AddModuleScore(Mets, features = list(dormancy_genes), name = "DormancyScore")
Mets <- AddModuleScore(Mets, features = list(stress_genes), name = "StressScore")

# Module scores will now be:
# ApoptosisScore1, DormancyScore1, StressScore1 in your metadata

# ---------------------------
# 3. Extract scores for IFNg-high and IFNg-low groups
# ---------------------------

# (assuming you already have IFNg_high_cells and IFNg_low_cells)

# Fetch scores
scores_low <- FetchData(Mets, vars = c("ProSurvivalScore1", "DormancyScore1", "StressScore1"), cells = IFNg_low_cells)
scores_high <- FetchData(Mets, vars = c("ProSurvivalScore1", "DormancyScore1", "StressScore1"), cells = IFNg_high_cells)

# Combine into one big dataframe
scores_df <- rbind(
  data.frame(scores_low, Group = "IFNg-low"),
  data.frame(scores_high, Group = "IFNg-high")
)

# ---------------------------
# 4. Plot all three programs
# ---------------------------

library(ggplot2)
library(reshape2)

# Melt the dataframe to long format for easier ggplot
scores_df_melt <- melt(scores_df, id.vars = "Group", variable.name = "Program", value.name = "Score")

# Violin plot for each score
ggplot(scores_df_melt, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white") +
  facet_wrap(~ Program, scales = "free_y") +
  theme_minimal() +
  ylab("Module Score") +
  xlab("") +
  ggtitle("Survival and Stress Programs in IFNg-high vs IFNg-low Cancer Cells")

# ----------------------------------------
# 1. Apoptosis score comparison (IFNg-low vs IFNg-high)
# ----------------------------------------
prosurvival_wilcox <- wilcox.test(scores_low$ProSurvivalScore1, scores_high$ProSurvivalScore1)
print(paste("Apoptosis Score p-value:", prosurvival_wilcox$p.value))

# ----------------------------------------
# 2. Dormancy score comparison (IFNg-low vs IFNg-high)
# ----------------------------------------
dormancy_wilcox <- wilcox.test(scores_low$DormancyScore1, scores_high$DormancyScore1)
print(paste("Dormancy Score p-value:", dormancy_wilcox$p.value))

# ----------------------------------------
# 3. Stress score comparison (IFNg-low vs IFNg-high)
# ----------------------------------------
stress_wilcox <- wilcox.test(scores_low$StressScore1, scores_high$StressScore1)
print(paste("Stress Score p-value:", stress_wilcox$p.value))

# ----------------------------------------
# 4. Autophagy score comparison (IFNg-low vs IFNg-high)
# ----------------------------------------
autophagy_wilcox <- wilcox.test(IFNlow_score$AutophagyScore1, IFNhigh_score$AutophagyScore1)
print(paste("Autophagy Score p-value:", autophagy_wilcox$p.value))

# Perform DEG analysis between IFNg-low and IFNg-high cells
deg_results <- FindMarkers(Mets, ident.1 = "IFNg-high", ident.2 = "IFNg-low", 
                           min.pct = 0.25, logfc.threshold = 0.25)

# Volcano plot of DEGs
write.csv(deg_results, file = "DEG_IFNgHigh_IFNgLow.csv")

deg_results$neg_adj.pVal <- -1*log(deg_results$p_val_adj, 10)
pThresh <- 10
plot(deg_results$avg_log2FC, deg_results$neg_adj.pVal, col=ifelse(deg_results$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 150), pch=ifelse(deg_results$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-10, 5, by=2.5), col="lightgray")
abline(h=seq(0,150, by=25), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(deg_results$avg_log2FC, deg_results$neg_adj.pVal, col=ifelse(deg_results$neg_adj.pVal > pThresh, ifelse(deg_results$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(deg_results$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(deg_results$neg_adj.pVal > pThresh, 1.4, 1.2))
ggplotly(ggplot(data = deg_results,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(deg_results), colour=ifelse(deg_results$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

#"Apoptosis Score p-value: 9.66694606247425e-23"
#"Dormancy Score p-value: 1.6887569953371e-59"
# "Stress Score p-value: 3.88830591033298e-12"
#"Autophagy Score p-value: 0.134034973915433"

# Assign a default label
Mets$IFNg_status <- "Mid"

# Label manually selected cells
Mets$IFNg_status[IFNg_low_cells] <- "IFNg-low"
Mets$IFNg_status[IFNg_high_cells] <- "IFNg-high"

# Set identities based on this column
Idents(Mets) <- "IFNg_status"

VlnPlot(Mets, features = c("Gadd45a", "Gadd45g", "Stat3", "Bcl2l11"))

# Set identities to IFNg status
Idents(Mets) <- "IFNg_status"

# Define the genes you're interested in
genes_of_interest <- c("Gadd45a", "Gadd45g", "Stat3", "Bcl2l11")

# Perform differential expression analysis
deg_results <- FindMarkers(Mets,
                           ident.1 = "IFNg-high",
                           ident.2 = "IFNg-low",
                           features = genes_of_interest)

# Print or view the p-values
deg_results[, c("avg_log2FC", "p_val", "p_val_adj")]


