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

Mets <- readRDS(file = "/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/RDS/Mets_only_post_clusterGESA.rds")

# ===== BC overlay on clusters =======

## ===== START of Lib.ID and LARRY.BC freq plotting =============================
### === Quantify Lib.BC-ID frequency across each HTO in Asc.met and Om.met 
### Subset the data
Idents(Mets) <- "Combined.HTO_group"

OmMet <- subset(Mets, idents = c("om:2_OmMet"))
AscMet <- subset(Mets, idents = c("asc:1_AscMet"))

# Extract metadata
OmMet.meta <- OmMet@meta.data
AscMet.meta <- AscMet@meta.data

##=== Function to create frequency table for each Combined.HTO_group
create_frequency_table <- function(meta_data, group_column, barcode_column) {
  meta_data %>%
    group_by(across(all_of(group_column)), across(all_of(barcode_column))) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = all_of(barcode_column), values_from = count, values_fill = list(count = 0))
}

# Create frequency tables
OmMet_freq_table <- create_frequency_table(OmMet.meta, "om_HTO_classification", "om_LibID_first_bc")
AscMet_freq_table <- create_frequency_table(AscMet.meta, "Asc_HTO_classification", "Asc_LibID_first_bc")

# Print the frequency tables
print(OmMet_freq_table)
print(AscMet_freq_table)

# plot heatmap of BC_libID
# Function to reshape and plot heatmap
plot_normalized_heatmap <- function(df, classification_col) {
  
  # Step 1: Remove not_detected
  df_clean <- df %>% select(-not_detected)
  
  # Step 2: Normalize per row (row-wise proportions)
  df_norm <- df_clean %>%
    column_to_rownames(classification_col) %>%
    as.matrix() %>%
    prop.table(1) %>%
    as_tibble(rownames = classification_col)
  
  # Step 3: Pivot for plotting
  df_long <- df_norm %>%
    pivot_longer(-all_of(classification_col), names_to = "Barcode", values_to = "Proportion")
  
  # Step 4: Sort Barcode by total frequency across all classifications
  desired_order <- c("R1_BC_ID6", "R1_BC_ID5", "R1_BC_ID4", "R1_BC_ID3", "R1_BC_ID2", "R1_BC_ID1")
  
  df_long <- df_long %>%
    filter(Barcode %in% desired_order) %>%
    mutate(Barcode = factor(Barcode, levels = desired_order))
  
  # Step 5: Plot heatmap
  ggplot(df_long, aes(x = Barcode, y = .data[[classification_col]], fill = Proportion)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "darkred", name = "Proportion") +
    theme_minimal() +
    labs(x = "Lib.BC-ID", y = "HTO Classification") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot heatmaps (normalized across raw HTO)
plot_normalized_heatmap(OmMet_freq_table, "om_HTO_classification")

pdf(file = "./Figures/Met_OmMet_freq_table_heatmap.pdf", width = 6, height = 2.3, paper='special')
plot_normalized_heatmap(OmMet_freq_table, "om_HTO_classification")
dev.off()

plot_normalized_heatmap(AscMet_freq_table, "Asc_HTO_classification")

pdf(file = "./Figures/Met_AscMet_freq_table_heatmap.pdf", width = 6, height = 3.5, paper='special')
plot_normalized_heatmap(AscMet_freq_table, "Asc_HTO_classification")
dev.off()

# Examine LARRY frequency
# Create frequency tables
OmMet_LARRY_freq_table <- create_frequency_table(OmMet.meta, "om_HTO_classification", "om_LARRY_first_bc")
AscMet_LARRY_freq_table <- create_frequency_table(AscMet.meta, "Asc_HTO_classification", "Asc_LARRY_first_bc")

OmMet_LARRY_freq_table
AscMet_LARRY_freq_table

# Plot heatmaps of Larry BC (normalized across raw HTO and sorted before plotting)
plot_normalized_heatmap_ggplot <- function(df, classification_col, top_n_barcodes = 50,
                                           exclude_barcodes = c("not_detected", "LARRY_ratio_less 2")) {
  
  # Step 1: Normalize per row
  df_norm <- df %>%
    column_to_rownames(classification_col) %>%
    as.matrix() %>%
    prop.table(1) %>%
    as_tibble(rownames = classification_col)
  
  # Step 2: Pivot long
  df_long <- df_norm %>%
    pivot_longer(-all_of(classification_col), names_to = "Barcode", values_to = "Proportion")
  
  # Step 3: Remove unwanted barcodes
  df_long <- df_long %>%
    filter(!Barcode %in% exclude_barcodes)
  
  # Step 4: Rank barcodes by total normalized frequency
  top_barcodes <- df_long %>%
    group_by(Barcode) %>%
    summarize(Total = sum(Proportion), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = top_n_barcodes) %>%
    pull(Barcode)
  
  # Step 5: Filter data for top barcodes and prepare for ggplot
  # Ensure the classification_col is a factor with original order if desired, or let ggplot handle alphabetically
  # Here we convert Barcode to a factor to maintain the desired order from top_barcodes
  df_plot <- df_long %>%
    filter(Barcode %in% top_barcodes) %>%
    mutate(Barcode = factor(Barcode, levels = top_barcodes)) %>%
    # Convert classification_col to a factor for desired row order in heatmap (if applicable)
    # By default, ggplot will arrange rows alphabetically. If a specific order is desired,
    # the classification_col also needs to be factored according to that order.
    # For now, we assume original order from the input df is retained for rows.
    mutate(!!sym(classification_col) := factor(!!sym(classification_col),
                                               levels = rev(unique(df_norm[[classification_col]])))) # Reverse to put first row at top
  
  # Step 6: Plot with ggplot2
  p <- ggplot(df_plot, aes(x = Barcode, y = !!sym(classification_col), fill = Proportion)) +
    geom_tile(color = "black", linewidth = 0.1) +
    # Use scale_fill_gradientn to specify colors at specific values
    # Colors will go from white (at 0) to a dark magenta (at max proportion)
    scale_fill_gradientn(colors = c("white", "darkmagenta"), # Start at white, end at darkmagenta
                         values = c(0, 1), # White at 0, darkmagenta at 1 (max proportion)
                         limits = c(0, max(df_plot$Proportion, na.rm = TRUE)), # Ensure scale starts at 0
                         name = "Normalized\nProportion",
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Normalized Proportion Heatmap:", classification_col),
         x = "Barcode",
         y = gsub("_", " ", classification_col)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  
  return(p) # ggplot objects can be returned and then saved
}

# To save a ggplot object to PDF:
plot_ommet_ggplot <- plot_normalized_heatmap_ggplot(OmMet_LARRY_freq_table, "om_HTO_classification", top_n_barcodes = 50)
plot_ommet_ggplot

ggsave(filename = "./Figures/Met_OmMet_LARRY_freq_table_heatmap_ggplot.pdf", 
       plot = plot_ommet_ggplot, width = 12, height = 1.65, units = "in")

plot_ascmet_ggplot <- plot_normalized_heatmap_ggplot(AscMet_LARRY_freq_table, "Asc_HTO_classification", top_n_barcodes = 50)
plot_ascmet_ggplot

ggsave(filename = "./Figures/Met_AscMet_LARRY_freq_table_heatmap_ggplot.pdf",
       plot = plot_ascmet_ggplot, width = 12, height = 2, units = "in")

## plot RAW barcode counts using based R + pheatmap ====
plot_raw_count_heatmap_ggplot <- function(df, classification_col, top_n_barcodes = 50,
                                          exclude_barcodes = c("not_detected", "LARRY_ratio_less 2")) {
  
  # Step 1: Prepare data with classification_col as rownames
  # For raw counts, we directly work with the input df (after ensuring it's numeric and clean)
  # Ensure the classification_col is not numeric columns that will be plotted as counts.
  df_counts <- df %>%
    # Make a copy to avoid modifying the original dataframe in the calling environment
    dplyr::select(-all_of(classification_col)) %>% # Exclude the classification column for rowname conversion
    as.matrix()
  
  # Ensure all columns are numeric
  if (!all(sapply(df_counts, is.numeric))) {
    stop("All columns (excluding the classification column) must be numeric for raw count heatmap.")
  }
  
  # Add classification_col back as a regular column for pivoting
  df_counts <- df %>%
    column_to_rownames(classification_col) %>%
    as.matrix() %>%
    as_tibble(rownames = classification_col)
  
  # Step 2: Pivot long
  df_long <- df_counts %>%
    pivot_longer(-all_of(classification_col), names_to = "Barcode", values_to = "Count")
  
  # Step 3: Remove unwanted barcodes
  df_long <- df_long %>%
    filter(!Barcode %in% exclude_barcodes)
  
  # Step 4: Rank barcodes by total raw count
  # We still rank by total sum, but now it's the sum of raw counts
  top_barcodes <- df_long %>%
    group_by(Barcode) %>%
    summarize(Total = sum(Count), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = top_n_barcodes) %>%
    pull(Barcode)
  
  # Step 5: Filter data for top barcodes and prepare for ggplot
  df_plot <- df_long %>%
    filter(Barcode %in% top_barcodes) %>%
    mutate(Barcode = factor(Barcode, levels = top_barcodes)) %>%
    mutate(!!sym(classification_col) := factor(!!sym(classification_col),
                                               levels = rev(unique(df_counts[[classification_col]])))) # Reverse to put first row at top
  
  # Step 6: Plot with ggplot2
  p <- ggplot(df_plot, aes(x = Barcode, y = !!sym(classification_col), fill = Count)) + # 'fill' is now 'Count'
    geom_tile(color = "black", linewidth = 0.1) + # Add white borders for cell separation
    # Modified color scale to ensure 0 is white
    scale_fill_gradientn(colors = c("white", "darkmagenta"), # Start at white, transition to darkmagenta
                         values = scales::rescale(c(0, max(df_plot$Count, na.rm = TRUE)), to = c(0, 1)), # Map 0 to 0 and max to 1
                         limits = c(0, max(df_plot$Count, na.rm = TRUE)), # Ensure scale starts at 0 and goes up to max observed count
                         name = "Raw\nCount", # Customize legend title for raw counts
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Raw Count Heatmap:", classification_col), # Update title
         x = "Barcode",
         y = gsub("_", " ", classification_col)) + # Clean up y-axis label
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8), # Rotate x-axis labels
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  return(p) # ggplot objects can be returned and then saved
}

plot_ommet_raw_ggplot <- plot_raw_count_heatmap_ggplot(OmMet_LARRY_freq_table, "om_HTO_classification", top_n_barcodes = 50)
plot_ommet_raw_ggplot

ggsave(filename = "./Figures/Met_OmMet_LARRY_freq_table_heatmap_raw_ggplot.pdf",
       plot = plot_ommet_raw_ggplot, width = 12, height = 1.65, units = "in")

plot_ascmet_raw_ggplot <- plot_raw_count_heatmap_ggplot(AscMet_LARRY_freq_table, "Asc_HTO_classification", top_n_barcodes = 50)
plot_ascmet_raw_ggplot

ggsave(filename = "./Figures/Met_AscMet_LARRY_freq_table_heatmap_raw_ggplot.pdf",
       plot = plot_ascmet_raw_ggplot, width = 12, height = 2, units = "in")

## ======== END of Lib.ID and LARRY.BC freq plotting ============================


# ===== Examine DEGs between Omental metastases and ascites ======= #############
Idents(Mets) <- "Combined.HTO_group"
DEG_Mets <- FindMarkers(Mets, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_Mets, file = "./Results/OmMet_compared_to_ascMet_DEG.csv")

DEG_Mets <- DEG_Mets %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj),
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 1  ~ "Up in Om",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Up in Asc",
      TRUE                               ~ "NS"
    )
  )

## Plot volcano
volcano_plot <- ggplot(DEG_Mets, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up in Om" = "red", "Up in Asc" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "Average log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "DEG Category"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# label top genes by effect size or p-value
top_genes <- DEG_Mets %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 2) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 30)  

volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = Inf
  )

print(volcano_plot)

pdf(file = "./Figures/Mets_DEG_asc-om_volcano.pdf", width = 8, height = 6, paper='special')
print(volcano_plot)
dev.off()

# ================ GESA of DEG Asc vs Om ===================================+===
# Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Step 1: Filter top significant genes for Asc and Om
deg_df <- DEG_Mets

# Define cutoffs
padj_cutoff <- 0.05
logfc_cutoff <- 0.5

# Up in Asc
up_in_asc <- deg_df %>%
  filter(significance == "Up in Asc", p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff)

# Up in Om
up_in_om <- deg_df %>%
  filter(significance == "Up in Om", p_val_adj < padj_cutoff, abs(avg_log2FC) > logfc_cutoff)

# Step 2: Map gene symbols to Entrez IDs
gene_to_entrez <- bitr(unique(c(up_in_asc %>% rownames(),
                                up_in_om %>% rownames())),
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Mm.eg.db)

# Merge Entrez IDs
up_in_asc_entrez <- gene_to_entrez %>%
  filter(SYMBOL %in% rownames(up_in_asc)) %>%
  pull(ENTREZID)

up_in_om_entrez <- gene_to_entrez %>%
  filter(SYMBOL %in% rownames(up_in_om)) %>%
  pull(ENTREZID)

# Step 3: Enrichment analysis using GO Biological Process
go_asc <- enrichGO(gene         = up_in_asc_entrez,
                   OrgDb        = org.Mm.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable     = TRUE)

go_om <- enrichGO(gene         = up_in_om_entrez,
                  OrgDb        = org.Mm.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

# Step 4: Visualize top enriched pathways
dotplot(go_asc, showCategory = 20, title = "GO BP Enrichment - Up in Asc")
dotplot(go_om, showCategory = 20, title = "GO BP Enrichment - Up in Om")

pdf(file = "./Figures/Mets_go_asc_enrichment.pdf", width = 8, height = 6, paper='special')
dotplot(go_asc, showCategory = 20, title = "GO BP Enrichment - Up in Asc")
dev.off()

pdf(file = "./Figures/Mets_go_om_enrichment.pdf", width = 8, height = 6, paper='special')
dotplot(go_om, showCategory = 20, title = "GO BP Enrichment - Up in Om")
dev.off()

# KEGG enrichment analysis
kegg_asc <- enrichKEGG(gene = up_in_asc_entrez,
                       organism = "mmu",
                       pvalueCutoff = 0.05)

kegg_om <- enrichKEGG(gene = up_in_om_entrez,
                      organism = "mmu",
                      pvalueCutoff = 0.05)

dotplot(kegg_asc, showCategory = 20, title = "KEGG - Up in Asc")
dotplot(kegg_om, showCategory = 20, title = "KEGG - Up in Om")

pdf(file = "./Figures/Mets_kegg_asc_enrichment.pdf", width = 8, height = 6, paper='special')
dotplot(kegg_asc, showCategory = 20, title = "KEGG BP Enrichment - Up in Asc")
dev.off()

pdf(file = "./Figures/Mets_kegg_om_enrichment.pdf", width = 8, height = 6, paper='special')
dotplot(kegg_om, showCategory = 20, title = "KEGG BP Enrichment - Up in Om")
dev.off()

## ===== GSEA enrichment of Asc and Om DEG =========
# STEP 1: Prepare ranking metric: signed -log10(padj)
deg_df <- DEG_Mets %>%
  mutate(SYMBOL = rownames(.),
         signed_stat = -log10(p_val_adj + 1e-300) * sign(avg_log2FC)) %>%
  arrange(desc(signed_stat))

# STEP 2: Convert SYMBOL to ENTREZID
symbol_to_entrez <- bitr(deg_df$SYMBOL,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Mm.eg.db)

# STEP 3: Join Entrez IDs and create geneList
deg_ranked <- deg_df %>%
  inner_join(symbol_to_entrez, by = "SYMBOL") %>%
  arrange(desc(signed_stat)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_list <- deg_ranked$signed_stat
names(gene_list) <- deg_ranked$ENTREZID

# Define input directory and list of gmt files
gmt_files <- list.files(path = "./assets/GSEA", pattern = "*.gmt", full.names = TRUE)

# Prepare ranked gene list
gene_list_symbol <- deg_df %>%
  filter(!is.na(signed_stat)) %>%
  arrange(desc(signed_stat)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_ranks <- gene_list_symbol$signed_stat
names(gene_ranks) <- gene_list_symbol$SYMBOL

# Loop through each gmt file
for (gmt_path in gmt_files) {
  gmt_name <- tools::file_path_sans_ext(basename(gmt_path))  # For naming output
  term2gene <- read.gmt(gmt_path)
  
  # Run GSEA
  gsea_result <- GSEA(geneList = gene_ranks,
                      TERM2GENE = term2gene,
                      pvalueCutoff = 0.01,
                      verbose = FALSE)
  
  # Annotate group based on NES
  gsea_result@result$group <- ifelse(gsea_result@result$NES > 0, "Om", "Asc")
  gsea_result@result$log10padj <- -log10(gsea_result@result$p.adjust + 1e-300)
  
  # Create dotplot
  p <- dotplot(gsea_result, showCategory = 10, split = "group") +
    facet_grid(. ~ group) +
    scale_color_gradient(
      low = "blue", high = "red",
      name = "-log10(adj p)",
      guide = "colorbar"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      strip.text.x = element_text(size = 12, face = "bold")
    )
  
  # Save to PDF
  output_file <- paste0("./Figures/GSEA_", gmt_name, "_dotplot_split_by_group.pdf")
  pdf(output_file, width = 7, height = 6)
  print(p)
  dev.off()
}

# === Examine IFNa and IFNg genesets ================
# Read hallmark gene set
hallmark_gmt <- read.gmt("./assets/GSEA/mh.all.v2024.1.Mm.symbols.gmt")

# Extract genes for the INTERFERON_ALPHA_RESPONSE set
ifna_genes <- hallmark_gmt %>%
  filter(term == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene)

# View
print(ifna_genes)

# Use your ifna_genes list and intersect with canonical IFN-gamma genes
canonical_ifng_genes <- c(
  "Irf1", "Stat1", "Cxcl9", "Cxcl10", "Cxcl11", "Tap1", "Tap2",
  "Psmb8", "Psmb9", "H2-K1", "H2-D1", "H2-Aa", "H2-Ab1", "H2-Eb1",
  "Cd74", "B2m", "Gbp2", "Gbp3", "Gbp5", "Nos2", "Ciita", "Socs1",
  "Casp1", "Casp8", "Ifi30", "Trim21", "Isg15"
)

# Overlap
overlap <- intersect(ifna_genes, canonical_ifng_genes)
length(overlap)
print(overlap)

Idents(Mets) <- "Combined.HTO_group"
DimPlot(Mets)
VlnPlot(Mets, features = overlap, assay = "SCT")

pdf(file = "./Figures/Mets_VlnPlot_IFNa-IFNg_overlap_Asc.VS.Om.pdf", width = 8, height = 10, paper='special')
VlnPlot(Mets, features = overlap, assay = "SCT")
dev.off()

Mets<- AddModuleScore(Mets, features = list(canonical_ifng_genes), name = "IFNG_score")
VlnPlot(Mets, features = "IFNG_score1", group.by = "Combined.HTO_group")

Mets<- AddModuleScore(Mets, features = list(ifna_genes), name = "IFNa_score")
VlnPlot(Mets, features = "IFNa_score1", group.by = "Combined.HTO_group")

VlnPlot(Mets, features = c("IFNa_score1", "IFNG_score1"), group.by = "Combined.HTO_group")

pdf(file = "./Figures/Mets_Asc-Om_Vln_IFN.score.pdf", width = 5, height = 6, paper='special')
VlnPlot(Mets, features = c("IFNa_score1", "IFNG_score1"), group.by = "Combined.HTO_group")
dev.off()


# === Examine the ractome ===
library(ReactomePA)

kegg_res <- gseKEGG(geneList = gene_list,
                    organism = "mmu",
                    pvalueCutoff = 0.05,
                    verbose = FALSE)

kegg_res@result$group <- ifelse(kegg_res@result$NES > 0, "Om", "Asc")

reactome_res <- gsePathway(geneList = gene_list,
                           organism = "mouse",
                           pvalueCutoff = 0.05,
                           verbose = FALSE)

reactome_res@result$group <- ifelse(reactome_res@result$NES > 0, "Om", "Asc")

dotplot(kegg_res, showCategory = 10, split = "group") +
  facet_grid(. ~ group) +
  scale_fill_gradient(low = "white", high = "#cc0033", name = "-log10(adj p)") +
  theme_bw(base_size = 12) +
  ggtitle("KEGG GSEA: Enriched Pathways in Asc vs Om")

pdf(file = "./Figures/Mets_Asc-Om_kegg_res.pdf", width = 8, height = 6, paper='special')
dotplot(kegg_res, showCategory = 10, split = "group") +
  facet_grid(. ~ group) +
  scale_fill_gradient(low = "white", high = "#cc0033", name = "-log10(adj p)") +
  theme_bw(base_size = 10) +
  ggtitle("KEGG GSEA: Enriched Pathways in Asc vs Om")
dev.off()

pdf(file = "./Figures/Mets_Asc-Om_reactome_res.pdf", width = 8, height = 6, paper='special')
dotplot(reactome_res, showCategory = 10, split = "group") +
  facet_grid(. ~ group) +
  scale_fill_gradient(low = "white", high = "#cc0033", name = "-log10(adj p)") +
  theme_bw(base_size = 12) +
  ggtitle("Reactome GSEA: Enriched Pathways in Asc vs Om")
dev.off()

# == Perform DEG analysis but with BC.ID filtering
# Compare highly enriched IDs vs low
# Check transcriptome differences among BC.IDs within individual treatment groups

Idents(Mets) <- "om_LibID_first_bc"
DEG_Om.Mets.BCID6.vs.BCID12 <- FindMarkers(Mets, ident.1 = c("R1_BC_ID6"), ident.2 = c("R1_BC_ID1", "R1_BC_ID2"), verbose = TRUE, logfc.threshold = 0.05)

Idents(Mets) <- "Asc_LibID_first_bc"
DEG_Asc.Mets.BCID6.vs.BCID12 <- FindMarkers(Mets, ident.1 = c("R1_BC_ID6"), ident.2 = c("R1_BC_ID1", "R1_BC_ID2"), verbose = TRUE, logfc.threshold = 0.05)

# no DEGs for Asc comparison (adjp is bad), but there are DEGs in Om - add to CRISPR list
write.csv(DEG_Om.Mets.BCID6.vs.BCID12, file = "./Results/DEG_Om.Mets.BCID6.vs.BCID12.csv")
write.csv(DEG_Asc.Mets.BCID6.vs.BCID12, file = "./Results/DEG_Asc.Mets.BCID6.vs.BCID12.csv")

# valcano plot DEG_Om.Mets.BCID6.vs.BCID12
DEG_Om.Mets.BCID6.vs.BCID12 <- DEG_Om.Mets.BCID6.vs.BCID12 %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj),
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 1  ~ "Up in BCID6",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Up in BCID1.2",
      TRUE                               ~ "NS"
    )
  )

# Plot volcano
volcano_plot <- ggplot(DEG_Om.Mets.BCID6.vs.BCID12, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Up in BCID6" = "red", "Up in BCID1.2" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "Average log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "DEG Category"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# Optional: label top genes by effect size or p-value
top_genes <- DEG_Om.Mets.BCID6.vs.BCID12 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 30)  

volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = Inf
  )

print(volcano_plot)

pdf(file = "./Figures/Mets_DEG_BCID1.2-BCID6_volcano.pdf", width = 8, height = 6, paper='special')
print(volcano_plot)
dev.off()

# === Plot LARRY BC freq ===
Idents(Mets) <- "Combined.HTO_group"
OmMet <- subset(Mets, idents = c("om:2_OmMet"))
OmMet@meta.data$BC_status <- ifelse(OmMet@meta.data$Combined.LARRY_group != "not_detected-om", "LARRY_detected", "LARRY_not_detected")
Idents(OmMet) <- "BC_status"
OmMet <- subset(OmMet, idents = c("LARRY_detected"))
Idents(OmMet) <- "Combined.LARRY_group"

# Draw a frequency map using ggplot
Om_id_freq <- table(OmMet$Combined.LARRY_group)
Om_freq_df <- data.frame(LARRY_BC = names(Om_id_freq), Frequency = as.integer(Om_id_freq))
Om_freq_df_sorted <- Om_freq_df[order(-Om_freq_df$Frequency), ]
Om_freq_df_sorted$LARRY_BC <- factor(Om_freq_df_sorted$LARRY_BC, levels = Om_freq_df_sorted$LARRY_BC[order(-Om_freq_df_sorted$Frequency)])
Om_freq_df_sorted <- Om_freq_df_sorted[-1, ]
# above code removes the first row which contains "LARRY_ratio_less 2 -om" frequencies instead of actual BCs
Om_freq_df_top50   <- Om_freq_df_sorted %>% 
  slice_max(Frequency, n = 50, with_ties = FALSE)

plot_Om_freq_df_top50 <- ggplot(Om_freq_df_top50, aes(x = LARRY_BC, y = Frequency)) + geom_bar(stat = "identity", fill = "brown3") + theme_minimal() + labs(title = "Frequency of Top50 LARRY BCs in OmMet", x = "LARRY_BC", y = "Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Om_freq_df_top50
ggsave(filename = "./Figures/Mets_OmMet_LARRY_freq_hist.pdf",
       plot = plot_Om_freq_df_top50, width = 8, height = 4, units = "in")

Idents(Mets) <- "Combined.HTO_group"
AscMet <- subset(Mets, idents = c("asc:1_AscMet"))
AscMet@meta.data$BC_status <- ifelse(AscMet@meta.data$Combined.LARRY_group != "not_detected-asc", "LARRY_detected", "LARRY_not_detected")
Idents(AscMet) <- "BC_status"
AscMet <- subset(AscMet, idents = c("LARRY_detected"))
Idents(AscMet) <- "Combined.LARRY_group"

# Draw a frequency map using ggplot
Asc_id_freq <- table(AscMet$Combined.LARRY_group)
Asc_freq_df <- data.frame(LARRY_BC = names(Asc_id_freq), Frequency = as.integer(Asc_id_freq))
Asc_freq_df_sorted <- Asc_freq_df[order(-Asc_freq_df$Frequency), ]
Asc_freq_df_sorted$LARRY_BC <- factor(Asc_freq_df_sorted$LARRY_BC, levels = Asc_freq_df_sorted$LARRY_BC[order(-Asc_freq_df_sorted$Frequency)])
Asc_freq_df_sorted <- Asc_freq_df_sorted[-1, ]
#above code removes the first row which contains "LARRY_ratio_less 2 -asc" frequencies instead of actual BCs
Asc_freq_df_top50   <- Asc_freq_df_sorted %>% 
  slice_max(Frequency, n = 50, with_ties = FALSE)

plot_Asc_freq_df_top50 <- ggplot(Asc_freq_df_top50, aes(x = LARRY_BC, y = Frequency)) + geom_bar(stat = "identity", fill = "brown3") + theme_minimal() + labs(title = "Frequency of Top50 LARRY BCs in AscMet", x = "LARRY_BC", y = "Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Asc_freq_df_top50
ggsave(filename = "./Figures/Mets_AscMet_LARRY_freq_hist.pdf",
       plot = plot_Asc_freq_df_top50, width = 8, height = 4, units = "in")


# ==== Define LARRY expansion clones  #####
# Bin top LARRY BCs and define expansion
Mets.meta <- Mets@meta.data
barcode_counts <- table(Mets.meta$Combined.LARRY_group)
Mets.meta$Expansion_value <- "Not_detected" # give a value first to all cells

expanded_barcodes <- names(barcode_counts[barcode_counts >= 6]) # use cut off or 6 for both Om and Asc
expanded_barcodes <- expanded_barcodes[startsWith(expanded_barcodes, "BC")]
expanded_barcodes # confirm

low_expanded_barcodes <- names(barcode_counts[barcode_counts >= 1 & barcode_counts < 6])
low_expanded_barcodes

not_detected_barcodes <- "LARRY_ratio_less 2-Mets"

# even the LARRY = 1, we considered sequencing and LARRY detection limit. 
# Thus, we still consider those  cells are "low expanded"
Mets.meta$Expansion_value[Mets.meta$Combined.LARRY_group %in% expanded_barcodes] <- "Expansion_high"
Mets.meta$Expansion_value[Mets.meta$Combined.LARRY_group %in% low_expanded_barcodes] <- "Expansion_low"
Mets.meta$Expansion_value[Mets.meta$Combined.LARRY_group %in% not_detected_barcodes] <- "Not_detected"

table(Mets.meta$Expansion_value)
colnames(Mets.meta)  

# Add the updated metadata back to the Seurat object
Mets <- AddMetaData(object = Mets, metadata = Mets.meta)

# Verify the new metadata slot
table(Mets@meta.data$Expansion_value)
Idents(Mets) <- "Expansion_value"

# Filter Seurat object to only include cells are not "Not_detected"
# Filter explicitly on both group and IFNa_score1 value
Mets_expansion_clean <- subset(Mets_expansion_filtered, subset = 
                                 Expansion_value %in% c("Expansion_high", "Expansion_low") &
                                 !is.na(IFNa_score1) & !is.na(IFNG_score1))
VlnPlot(Mets_expansion_clean, features = c("IFNa_score1", "IFNG_score1"), group.by = "Expansion_value")

# plot IFN score
pdf(file = "./Figures/Mets_IFNscore_expansion_high-low.pdf", width = 6, height = 5, paper='special')
VlnPlot(Mets_expansion_clean, features = c("IFNa_score1", "IFNG_score1"), group.by = "Expansion_value")
dev.off()

# test overall P value between groups. 
# Define features to test
features_to_test <- c("IFNa_score1", "IFNG_score1")

# Initialize results list
pval_results <- list()

# Loop through features
for (feature in features_to_test) {
  data_values <- FetchData(Mets_expansion_clean, vars = c(feature, "Expansion_value"))
  
  test_result <- wilcox.test(
    data_values[[feature]] ~ data_values$Expansion_value,
    exact = FALSE
  )
  
  # Store results
  pval_results[[feature]] <- data.frame(
    Feature = feature,
    p_value = test_result$p.value,
    W_statistic = test_result$statistic
  )
}

# Combine and save
pval_df <- bind_rows(pval_results)
write.csv(pval_df, file = "./Results/IFNscore_Expansion_WilcoxP.csv", row.names = FALSE)

# save RDS
saveRDS(Mets, file = "/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/RDS/Mets_only_post_clusterGESA_with_expansion.rds")

# == Mets <- readRDS(file = "/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/RDS/Mets_only_post_clusterGESA_with_expansion.rds")

Idents(Mets) <- "Combined.HTO_group"

OmMet <- subset(Mets, idents = c("om:2_OmMet"))
AscMet <- subset(Mets, idents = c("asc:1_AscMet"))

OmMet_filtered <- subset(Mets, subset = om_LibID_first_bc != "not_detected")
VlnPlot(OmMet_filtered, features = c("IFNa_score1", "IFNG_score1"), group.by = "om_LibID_first_bc") 

pdf(file = "./Figures/OmMet_IFNscore_LibID1-6.pdf", width = 6, height = 5)
VlnPlot(OmMet_filtered, features = c("IFNa_score1", "IFNG_score1"), group.by = "om_LibID_first_bc") 
dev.off()

AscMet_filtered <- subset(AscMet, subset = Asc_LibID_first_bc != "not_detected")
VlnPlot(AscMet_filtered, features = c("IFNa_score1", "IFNG_score1"), group.by = "Asc_LibID_first_bc") 

pdf(file = "./Figures/AscMet_IFNscore_LibID1-6.pdf", width = 6, height = 5)
VlnPlot(AscMet_filtered, features = c("IFNa_score1", "IFNG_score1"), group.by = "Asc_LibID_first_bc") 
dev.off()


# == DEG OmMet object. NOTE before bigger constrast, when doing DEG, 
# we define high.expansion >=6 and low.expansion as BC < 3
# Bin top LARRY BCs and define expansion
OmMet.meta <- OmMet@meta.data
barcode_counts <- table(OmMet.meta$Combined.LARRY_group)
OmMet.meta$Expansion_value <- "Not_detected" # give a value first to all cells

expanded_barcodes <- names(barcode_counts[barcode_counts >= 6]) # use cut off or 6 for both Om and Asc
expanded_barcodes <- expanded_barcodes[startsWith(expanded_barcodes, "BC")]
expanded_barcodes # confirm

low_expanded_barcodes <- names(barcode_counts[barcode_counts >= 1 & barcode_counts < 3])
low_expanded_barcodes

not_detected_barcodes <- "LARRY_ratio_less 2-OmMet"

# even the LARRY = 1, we considered sequencing and LARRY detection limit. 
# Thus, we still consider those  cells are "low expanded"
OmMet.meta$Expansion_value[OmMet.meta$Combined.LARRY_group %in% expanded_barcodes] <- "Expansion_high"
OmMet.meta$Expansion_value[OmMet.meta$Combined.LARRY_group %in% low_expanded_barcodes] <- "Expansion_low"
OmMet.meta$Expansion_value[OmMet.meta$Combined.LARRY_group %in% not_detected_barcodes] <- "Not_detected"

table(OmMet.meta$Expansion_value)
colnames(OmMet.meta)

# Add the updated metadata back to the Seurat object
OmMet <- AddMetaData(object = OmMet, metadata = OmMet.meta)

# Verify the new metadata slot
table(OmMet@meta.data$Expansion_value)
Idents(OmMet) <- "Expansion_value"

# Repeat with AscMet object 
# Bin top LARRY BCs and define expansion
AscMet.meta <- AscMet@meta.data
barcode_counts <- table(AscMet.meta$Combined.LARRY_group)
AscMet.meta$Expansion_value <- "Not_detected" # give a value first to all cells

expanded_barcodes <- names(barcode_counts[barcode_counts >= 6]) # use cut off or 6 for both Om and Asc
expanded_barcodes <- expanded_barcodes[startsWith(expanded_barcodes, "BC")]
expanded_barcodes # confirm

low_expanded_barcodes <- names(barcode_counts[barcode_counts >= 1 & barcode_counts < 3])
low_expanded_barcodes

not_detected_barcodes <- "LARRY_ratio_less 2-AscMet"

# even the LARRY = 1, we considered sequencing and LARRY detection limit. 
# Thus, we still consider those  cells are "low expanded"
AscMet.meta$Expansion_value[AscMet.meta$Combined.LARRY_group %in% expanded_barcodes] <- "Expansion_high"
AscMet.meta$Expansion_value[AscMet.meta$Combined.LARRY_group %in% low_expanded_barcodes] <- "Expansion_low"
AscMet.meta$Expansion_value[AscMet.meta$Combined.LARRY_group %in% not_detected_barcodes] <- "Not_detected"

table(AscMet.meta$Expansion_value)
colnames(AscMet.meta)

# Add the updated metadata back to the Seurat object
AscMet <- AddMetaData(object = AscMet, metadata = AscMet.meta)

# Verify the new metadata slot
table(AscMet@meta.data$Expansion_value)
Idents(AscMet) <- "Expansion_value"


# === DEGs between expanded and non-expanded clones
DEG_Om.Met_expansion <- FindMarkers(OmMet, ident.1 = "Expansion_high", ident.2 = "Expansion_low", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_Om.Met_expansion, file = "./Results/DEG_Om.Met_expansion_high-vs-low.csv")

# Filter for significant upregulated genes in both expansion and BCID in OmMet
sig_up_OmMet_expansion <- DEG_Om.Met_expansion[DEG_Om.Met_expansion$p_val_adj < 0.05 & DEG_Om.Met_expansion$avg_log2FC > 0, ]
sig_up_OmMet_BCID6_vs_BCID12 <- DEG_Om.Mets.BCID6.vs.BCID12[DEG_Om.Mets.BCID6.vs.BCID12$p_val_adj < 0.05 & DEG_Om.Mets.BCID6.vs.BCID12$avg_log2FC > 0, ]

# Find the common upregulated genes
common_genes <- intersect(rownames(sig_up_OmMet_expansion), rownames(sig_up_OmMet_BCID6_vs_BCID12))

# Extract the common genes from the first dataset (you can choose either)
common_genes_data <- sig_up_OmMet_expansion[common_genes, ]

# Write to CSV
write.csv(common_genes_data, "./Results/OmMet_common_upregulated_genes.csv")

# ==== compare DEGs between cells of Mets@meta.data$om_LibID_first_bc == R1_BC_ID6 and Mets@meta.data$Asc_LibID_first_bc == R1_BC_ID6  cells. 
# Create a new identity class in the Seurat object to define these two populations
Mets$LibID_bc_group <- NA  # initialize new metadata column

# Assign group labels
Mets$LibID_bc_group[Mets@meta.data$om_LibID_first_bc == "R1_BC_ID6"] <- "Om_R1_BC_ID6"
Mets$LibID_bc_group[Mets@meta.data$Asc_LibID_first_bc == "R1_BC_ID6"] <- "Asc_R1_BC_ID6"

# Confirm group sizes
table(Mets$LibID_bc_group)
# Should return something like:
# Asc_R1_BC_ID6  Om_R1_BC_ID6 
#            39           268 

# Set identities based on this grouping
Idents(Mets) <- "LibID_bc_group"

# Filter Seurat object to only include cells with non-NA group assignments
Mets_filtered <- subset(Mets, subset = !is.na(LibID_bc_group))

# Set identities again
Idents(Mets_filtered) <- "LibID_bc_group"

# Plot violin plot
VlnPlot(Mets_filtered, features = c("IFNa_score1", "IFNG_score1"), group.by = "LibID_bc_group")

pdf(file = "./Figures/Mets_IFNscore_LibID6_byOmAsc_filtered.pdf", width = 6, height = 5)
VlnPlot(Mets_filtered, features = c("IFNa_score1", "IFNG_score1"), group.by = "LibID_bc_group")
dev.off()

# test overall P value between groups. 

# Define features to test
features_to_test <- c("IFNa_score1", "IFNG_score1")

# Initialize results list
pval_results_LibID_bc <- list()

# Loop through features
for (feature in features_to_test) {
  data_values <- FetchData(Mets_filtered, vars = c(feature, "LibID_bc_group"))
  
  test_result <- wilcox.test(
    data_values[[feature]] ~ data_values$LibID_bc_group,
    exact = FALSE
  )
  
  # Store results
  pval_results_LibID_bc[[feature]] <- data.frame(
    Feature = feature,
    p_value = test_result$p.value,
    W_statistic = test_result$statistic
  )
}

# Combine and save
pval_LibID_bc_df <- bind_rows(pval_results_LibID_bc)
write.csv(pval_LibID_bc_df, file = "./Results/IFNscore_LibID_WilcoxP.csv", row.names = FALSE)


# Perform DEG analysis
DEG_LibID6_Om.Asc <- FindMarkers(Mets,
                                 ident.1 = "Om_R1_BC_ID6",
                                 ident.2 = "Asc_R1_BC_ID6",
                                 assay = "SCT",    
                                 logfc.threshold = 0.5,
                                 min.pct = 0.1)

# View top DEGs
head(DEG_LibID6_Om.Asc)

## Plot volcano
DEG_LibID6_Om.Asc <- DEG_LibID6_Om.Asc %>%
  mutate(
    p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 1  ~ "Up in Om.BCID6",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Up in Asc.BCID6",
      TRUE                               ~ "NS"
    )
  )
DEG_LibID6_Om.Asc$gene <- rownames(DEG_LibID6_Om.Asc)
write.csv(DEG_LibID6_Om.Asc, file = "./Results/DEG_LibID6_Om-Asc.csv")


volcano_plot <- ggplot(DEG_LibID6_Om.Asc, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up in Om.BCID6" = "red", "Up in Asc.BCID6" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot DEG_LibID6_Om.Asc",
    x = "Average log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "DEG Category"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# label top genes by effect size or p-value
top_genes <- DEG_LibID6_Om.Asc %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 2) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 30)  

# Define genes to label
genes_to_label <- c("Gbp2b", "Marco", "Ifngr1", "Ms4a8a")

# Create label data frames
top_genes <- DEG_LibID6_Om.Asc %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 2) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 30)

manual_genes <- DEG_LibID6_Om.Asc %>%
  filter(gene %in% genes_to_label)

# Combine and remove duplicates based on gene name
genes_to_highlight <- bind_rows(top_genes, manual_genes) %>%
  distinct(gene, .keep_all = TRUE)

volcano_plot <- volcano_plot +
  geom_text_repel(
    data = genes_to_highlight,
    aes(label = rownames(genes_to_highlight)),
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = Inf
  )

print(volcano_plot)

pdf(file = "./Figures/Mets_DEG_LibID6_Om-Asc_volcano.pdf", width = 6, height = 5, paper='special')
print(volcano_plot)
dev.off()


# Load the library
library(ggvenn)

# Use the 'gene_list' created in the previous method.
# If you skipped it, here is the setup again:
genes_expansion <- rownames(sig_up_OmMet_expansion)
genes_BCID_comparison <- rownames(sig_up_OmMet_BCID6_vs_BCID12)

gene_list <- list(
  `Expansion` = genes_expansion,
  `BCID6 vs BCID12` = genes_BCID_comparison
)

# Create the plot
# ggvenn works directly with the named list.
p <- ggvenn(
  gene_list, 
  fill_color = c("#ff9966", "#21908dff"), # Set fill colors
  stroke_size = 0.5,                       # Circle border size
  set_name_size = 4,                       # Text size for set names
  text_size = 5                            # Text size for counts
)

# Display the plot
print(p)

ggsave("./Figures/Mets_Expansion_Venn.pdf", plot = p, width = 8, height = 6)


#=== overlay Expanded LARRY BC on cluster 
Idents(Mets) <- "seurat_clusters"
DimPlot(Mets, split.by = "Expansion_value", label = TRUE)

pdf(file = "./Figures/Mets_DimPlot_expansion_lable.pdf", width = 10, height = 4, paper='special')
DimPlot(Mets, split.by = "Expansion_value")
dev.off()

# Subset the Seurat object to include only the cells in the 3 specified groups
# subset_cells <- subset(Mets, Combined.LARRY_group %in% c("BC-1-om", "BC-4-om", "BC-49-om",  "BC-75-om", "BC-330-om", "BC-483-om"))
# only subset top 3 the clarity  
subset_cells <- subset(Mets, Combined.LARRY_group %in% 
                         c("BC-1-om", "BC-4-om", "BC-49-om"))

# Set the identities of the cells to the 'Combined.LARRY_group'
Idents(subset_cells) <- "Combined.LARRY_group"

# Plot the UMAP with different colors for the three cell populations
DimPlot(subset_cells, reduction = "umap", group.by = "Combined.LARRY_group", 
        cols = c("BC-1-om" = "red", "BC-4-om" = "blue", "BC-49-om" = "#ffcc33", 
                   "BC-75-om" = "green", "BC-330-om" = "purple", "BC-483-om" = "orange"),
        pt.size = 3) +
  labs(title = "Top expanded clones")

DimPlot(Mets, label = TRUE)

# Extract cluster identities
cluster_ids <- Idents(Mets)

# Manually define Seurat cluster colors (based on the default DimPlot colors)
seurat_cluster_colors <- c(
  "0" = "#F8766D",  # Adjust based on Seurat default colors
  "1" = "#7CAE00",
  "2" = "#619CFF",
  "3" = "#C77CFF",
  "4" = "#E76BF3",
  "5" = "#FF8C00"
)

# plot with same xlim ylim 
# Extract UMAP embeddings from both objects
umap_subset <- Embeddings(subset_cells, reduction = "umap")
umap_Mets <- Embeddings(Mets, reduction = "umap")

# Determine shared axis limits
x_limits <- range(c(umap_subset[, 1], umap_Mets[, 1]))
y_limits <- range(c(umap_subset[, 2], umap_Mets[, 2]))

# First UMAP plot
p1 <- DimPlot(subset_cells, reduction = "umap", group.by = "Combined.LARRY_group", 
              cols = c("BC-1-om" = "red", "BC-4-om" = "blue", "BC-49-om" = "#ffcc33", 
                       "BC-75-om" = "green", "BC-330-om" = "purple", "BC-483-om" = "orange"),
              pt.size = 2) +
  xlim(x_limits) + ylim(y_limits) +
  labs(title = "Top expanded clones")

# Second UMAP plot
p2 <- DimPlot(Mets, label = FALSE, pt.size = 2) +
  xlim(x_limits) + ylim(y_limits) +
  labs(title = "Mets Cells")

# Display both plots side-by-side
p1 + p2

# Save the combined plot (p1+p2) as a PDF
ggsave("./Figures/Mets_combined_umap_with_top_LARRY.pdf", plot = p1 + p2, width = 10, height = 4.5, device = "pdf")

hallmark_features_to_plot <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE11",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE12",
  "HALLMARK_MITOTIC_SPINDLE1",
  "HALLMARK_HYPOXIA6"
)

p_expand <- FeaturePlot(
  object     = subset_cells,
  features   = hallmark_features_to_plot,
  ncol       = 2,
  cols       = c("white", "#cc6633"),
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size    = 0.7,      # how big each dot is
  order      = TRUE      # draw high‐expressing cells on top
)

print(p_expand)

pdf(file = "./Figures/Mets_FeaturePlot_top_hallmark_features_expanded_LARRY.pdf", width = 8.5, height = 8, paper='special')
print(p_expand)
dev.off()

Idents(Mets) <- "Combined.HTO_group"
DimPlot(Mets)
VlnPlot(Mets, features = c("Mki67", "Top2a"))
pdf(file = "./Figures/Mets_FeaturePlot_Mki67-Top2a.pdf", width = 5, height = 4, paper='special')
VlnPlot(Mets, features = c("Mki67", "Top2a"))
dev.off()

pdf(file = "./Figures/Mets_FeaturePlot_Gbp2b-Marco-Ifngr1.pdf", width = 7, height = 4, paper='special')
VlnPlot(Mets, features = c("Gbp2b", "Marco","Ifngr1"))
dev.off()

FeaturePlot(Mets, features = c("Gbp2b", "Marco","Ifngr1", "Ms4a8a"), order = TRUE)

## === Frequency plot for Combined.LARRY_group on Seurat cluster ==#############
# Create a frequency table for Combined.LARRY_group and seurat_clusters
freq_table <- table(subset_cells$Combined.LARRY_group, subset_cells$seurat_clusters)

# Convert the table to a data frame for easy plotting
freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("Population", "SeuratCluster", "Frequency")

# Create a frequency table for Combined.LARRY_group and seurat_clusters
freq_table <- table(subset_cells$Combined.LARRY_group, subset_cells$seurat_clusters)

# Generate UMAP plot with fixed colors
Idents(Mets) <- "seurat_clusters"
p3 <- DimPlot(Mets, label = TRUE, cols = seurat_cluster_colors) +
  labs(title = "UMAP of Seurat Clusters")

seurat_cluster_colors <- c(
  "0" = "#F8766D",  # Adjust based on Seurat default colors
  "1" = "#7CAE00",
  "2" = "#619CFF",
  "3" = "#C77CFF",
  "4" = "#E76BF3",
  "5" = "#FF8C00"
)

# Ensure Population is sorted by Frequency in descending order
freq_df <- freq_df %>%
  arrange(desc(Frequency)) %>%
  mutate(Population = factor(Population, levels = unique(Population)))  # Maintain sorted order
p3

# Generate the bar plot with sorted Population order, rotated x-axis labels, and a box around the plot
p4 <- ggplot(freq_df, aes(x = Population, y = Frequency, fill = SeuratCluster)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(title = "Frequency of Seurat Clusters by Population", 
       x = "Population", 
       y = "Frequency") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add black box around the plot
  ) +
  scale_fill_manual(values = seurat_cluster_colors)  # Ensure color consistency with UMAP
p4

# Display both plots side by side
p3 + p4

# Save both plots in a single PDF file
ggsave("./Figures/Mets_LARRY_combined_umap_freq_plots_p3-4.pdf", plot = p3 + p4, width = 8, height = 4, device = "pdf")

# ==== Check CRISPR hits via violin plot + FOLR1 expression on cancer cells ====
Mets@meta.data$Combined.HTO_group.order <- recode(Mets@meta.data$Combined.HTO_group, "asc:1_AscMet" = "01_AscMet", "om:2_OmMet" = "02_OmMet")
Idents(Mets) <- "Combined.HTO_group.order"
VlnPlot(Mets, features = c("Ccr1", "Cd274"))

# Specify the order of the groups for plotting
group_order <- c("01_AscMet", "02_OmMet")

# Plot the VlnPlot with the specified group order and gene expression for Ifng
VlnPlot(Mets, features = c("Gbp2b", "Marco", "Ms4a8a", "Cd300lb"), group.by = "Combined.HTO_group.order", 
        idents = group_order, ncol = 2, 
        pt.size = 1)

# Subset the Seurat object to only include AscMet and OmMet
Mets_subset <- subset(Mets, Combined.HTO_group.order %in% c("01_AscMet", "02_OmMet"))

# Get the number of cells in the smaller group (AscMet)
n_cells_AscMet <- sum(Mets_subset$Combined.HTO_group.order == "01_AscMet")

# Get the cells from the OmMet group and downsample to match AscMet
OmMet_cells <- WhichCells(Mets_subset, ident = "02_OmMet")

# Randomly sample the same number of cells as AscMet
set.seed(42)  # Set a seed for reproducibility
OmMet_cells_downsampled <- sample(OmMet_cells, size = n_cells_AscMet)

# Subset the Seurat object to only include the downsampled OmMet cells and all AscMet cells
Mets_balanced <- subset(Mets_subset, cells = c(WhichCells(Mets_subset, ident = "01_AscMet"), OmMet_cells_downsampled))

# Now, plot the VlnPlot for Ifng expression
VlnPlot(Mets_balanced, features =  c("Gbp2b", "Marco", "Ms4a8a", "Cd300lb"), group.by = "Combined.HTO_group.order", 
        pt.size = 1)

VlnPlot(Mets_balanced, features =  c("Ifngr1"), group.by = "Combined.HTO_group.order", pt.size = 1)

## ====  Plot specific gene with p-value=====
# Extract normalized expression matrix from SCT assay
gene_expression_data <- Mets_balanced[["SCT"]]@data

# Gene of interest
gene_of_interest <- "Gbp2b"

# Check if gene exists
if (gene_of_interest %in% rownames(gene_expression_data)) {
  
  # Get expression values for the gene
  expr_values <- gene_expression_data[gene_of_interest, ]
  
  # Get metadata and extract group information
  group_info <- Mets_balanced@meta.data$Combined.HTO_group.order
  names(group_info) <- colnames(Mets_balanced)
  
  # Create a data frame with expression and group labels
  df <- data.frame(
    expression = as.numeric(expr_values),
    group = as.factor(group_info)
  )
  
  # Subset data for the two groups to compare
  df_subset <- df[df$group %in% c("01_AscMet", "02_OmMet"), ]
  
  # Perform Wilcoxon rank-sum test
  wilcox_result <- wilcox.test(expression ~ group, data = df_subset)
  
  # Print test result
  print(wilcox_result)
  
} else {
  message("Gene not found in the dataset.")
}

p_val <- signif(wilcox_result$p.value, 3)
label_text <- paste0("p = ", p_val)

ggplot(df_subset, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = paste("Wilcoxon test for", gene_of_interest)) +
  annotate("text", x = 1.5, y = max(df_subset$expression, na.rm = TRUE) * 1.05, 
           label = label_text, size = 5)

Idents(Mets) <- "seurat_clusters"
DimPlot(Mets)

saveRDS(Mets, file = "/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/RDS/Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")

Mets <- readRDS(file = "/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/RDS/Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")

# ==== Mann-Whitney-Wilcoxon Gene Set Test (MWW-GST) analysis === #
devtools::install_github("YosefLab/VISION")
devtools::install_github("BorchLab/escape")

library(escape)
# Fetch mouse Hallmark gene sets
gs.hallmark.mouse <- getGeneSets(
  species = "Mus musculus",
  library = "H"
)

# compute per-cell UCell scores and store them as a new assay
Mets <- runEscape(
  input.data     = Mets,
  gene.sets      = gs.hallmark.mouse,
  method         = "UCell",
  min.size       = 5,
  new.assay.name = "escape.UCell_H",
  sets.size      = 2000,
  sets.n         = 5
)

nes_matrix <- GetAssayData(Mets, assay = "escape.UCell_H", slot = "data")
dim(nes_matrix)          # (#signatures × #cells)
head(nes_matrix[, 1:5])  # first few cells

# Visualize pathway activity 
DefaultAssay(Mets) <- "escape.UCell_H"
FeaturePlot(
  Mets,
  features  = "HALLMARK-INTERFERON-ALPHA-RESPONSE",
  reduction = "umap"
)

FeaturePlot(
  Mets,
  features  = "HALLMARK-INTERFERON-GAMMA-RESPONSE",
  reduction = "umap"
)

h <- Mets
h <- ScaleData(h, features = rownames(nes_matrix))
h <- RunPCA(h, features = rownames(nes_matrix), verbose = FALSE)
h <- FindNeighbors(h, dims = 1:10)
h <- FindClusters(h, resolution = 0.5)
h <- RunUMAP(h, dims = 1:10)
DimPlot(h, label = TRUE, split.by = "Combined.HTO_group") + ggtitle("Clusters by Hallmark NES")

pdf(file = "./Figures/Met_escape_enrichment_clusters.pdf", width = 10, height = 6, paper='special')
DimPlot(h, label = TRUE, split.by = "Combined.HTO_group") + ggtitle("Clusters by Hallmark NES")
dev.off()

# This will test for "marker pathways" per cluster
DefaultAssay(h) <- "escape.UCell_H"
hallmark_markers <- FindAllMarkers(
  h,
  only.pos     = TRUE,      # only positive markers (enriched)
  assay        = "escape.UCell_H",
  slot         = "data",    # use normalized NES scores
  min.pct      = 0.1,
  logfc.threshold = 0.1,
  test.use     = "wilcox"   # nonparametric (default is fine for NES)
)

# See the top markers per cluster:
top_hallmarks_clustermarker <- hallmark_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(5, avg_log2FC)

VlnPlot(
  h,
  features = "HALLMARK-INTERFERON-ALPHA-RESPONSE",  # change to pathway of interest
  group.by = "seurat_clusters",
  assay    = "escape.UCell_H",
  slot     = "data",
  pt.size  = 0.5  # no points, for cleaner plot
) + ggtitle("IFN-a Response NES across clusters")

pdf(file = "./Figures/Met_escape_enrichment_HALLMARK-INTERFERON-ALPHA-RESPONSE.pdf", width = 6, height = 6, paper='special')
VlnPlot(
  h,
  features = "HALLMARK-INTERFERON-ALPHA-RESPONSE",  # change to pathway of interest
  group.by = "seurat_clusters",
  assay    = "escape.UCell_H",
  slot     = "data",
  pt.size  = 0.5  # no points, for cleaner plot
) + ggtitle("IFN-a Response NES across clusters")
dev.off()

top_pathways <- hallmark_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(3, avg_log2FC) %>%
  dplyr::pull(gene) %>%
  unique()

VlnPlot(
  h,
  features = top_pathways,
  group.by = "seurat_clusters",
  assay    = "escape.UCell_H",
  slot     = "data",
  pt.size  = 0.2,
  stack    = FALSE,
  flip     = FALSE
)

pdf(file = "./Figures/Met_escape_enrichment_escape.UCell_H.pdf", width = 16, height = 20, paper='special')
VlnPlot(
  h,
  features = top_pathways,
  group.by = "seurat_clusters",
  assay    = "escape.UCell_H",
  slot     = "data",
  pt.size  = 0.2,
  stack    = FALSE,
  flip     = TRUE
)
dev.off()

pdf(file = "./Figures/Met_escape_enrichment_escape.UCell_H_stacked.pdf", width = 10, height = 20, paper='special')
VlnPlot(
  h,
  features = top_pathways,
  group.by = "seurat_clusters",
  assay    = "escape.UCell_H",
  slot     = "data",
  pt.size  = 0.2,
  stack    = TRUE,
  flip     = TRUE
)
dev.off()

pdf(file = "./Figures/Met_escape_enrichment_escape.UCell_heatmap.pdf", width = 8, height = 5, paper='special')
DoHeatmap(
  h,
  features = top_pathways,
  assay    = "escape.UCell_H",
  slot     = "data",
  group.by = "seurat_clusters"
) +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    limits = c(-0.3, 0.3) # adjust if needed based on your NES range
  ) +
  ggtitle("Hallmark NES Heatmap")
dev.off()

# === performing CytoTRACE analysis ===
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") #installing
library(CytoTRACE2) 

# For Seurat v5 and later:
counts_matrix <- GetAssayData(Mets, assay = "RNA", layer = "counts")

# Run CytoTRACE2
ct2_results <- cytotrace2(counts_matrix)

# The results are stored in ct2_results$score, which is named by cell
# Make sure the rownames of ct2_results match the cell names in Mets
head(rownames(ct2_results))   # Should look like "Om_AAACCCAAGCGTCTCG-1", etc.
head(colnames(Mets))          # Should be same format

# Create a vector of CytoTRACE2 scores named by cell barcode
cytotrace_scores <- ct2_results$CytoTRACE2_Score
names(cytotrace_scores) <- rownames(ct2_results)

# Assign CytoTRACE2 scores to the Seurat object's metadata
Mets$CytoTRACE2 <- cytotrace_scores[colnames(Mets)]
Mets$CytoTRACE2_Potency <- ct2_results$CytoTRACE2_Potency[colnames(Mets)]
Mets$CytoTRACE2_Relative <- ct2_results$CytoTRACE2_Relative[colnames(Mets)]

# UMAP plot colored by CytoTRACE2 scores
pdf(file = "./Figures/Met_CytoTRACE2_featureplot.pdf", width = 5, height = 5, paper='special')
FeaturePlot(Mets, features = "CytoTRACE2", cols = c("white", "red"))
dev.off()

# Violin plot by cluster
pdf(file = "./Figures/Met_CytoTRACE2_vlnplot.pdf", width = 6, height = 5, paper='special')
VlnPlot(Mets, features = "CytoTRACE2", group.by = "seurat_clusters")
dev.off()

Idents(Mets) <- "seurat_clusters"
DimPlot(Mets)

colnames((Mets@meta.data))

# == save Mets dataset as Mets.h5ad for Cellrank2

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)

object.i <- Mets
object.i[["RNA"]] <- CreateAssayObject(counts = object.i[["RNA"]]$counts) #Convert Assay v5 to Assay
SeuratDisk::SaveH5Seurat(object.i, filename = "/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/Looms/Mets.h5Seurat", overwrite = TRUE)
SeuratDisk::Convert("/Users/S208205/Library/CloudStorage/Box-Box/S.Zhang_Lab/UTSW_DATA/R02_EA/MetTag_BC_OVA/DATA/Looms/Mets.h5Seurat", dest = "h5ad", overwrite = TRUE)


