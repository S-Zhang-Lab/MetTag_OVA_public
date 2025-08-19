#### SF METHOD for GSEA #

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
setwd("C:/Users/Zhang Lab/Desktop/Emma/R01_and_2_w_LARRY")

# Load required packages specific for pathway enrichment
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)

# DEG analysis between OmMet and AscMet
Mets <- readRDS("Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")
Idents(Mets) <- "Combined.HTO_group"
DEG_OmMet_vs_AscMet <- FindMarkers(
  object = Mets,
  ident.1 = "om:2_OmMet",
  ident.2 = "asc:1_AscMet",
  verbose = TRUE,
  min.pct = 0.25,            # Minimum percentage of cells expressing the gene
  logfc.threshold = 0.25     # Log fold-change threshold
)

# View the top DEGs 
head(DEG_OmMet_vs_AscMet)
head(DEG_Trem2)

### Repeat the same procedure for Gbp high vs Gbp low cells
DEG_Gbp2b <- read.csv("DEG_Gbp2bHigh_vs_Low.csv", header = TRUE, stringsAsFactors = FALSE)

### GSEA Preparation #####
# Sort by log fold change to get up-regulated and down-regulated genes
upregulated_genes <- DEG_OmMet_vs_AscMet %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(500)

downregulated_genes <- DEG_OmMet_vs_AscMet %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  head(500)

upregulated_genes_AM <- DEG_Macs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%  # Added p-value filtering
  arrange(desc(avg_log2FC)) %>%
  head(300)

upregulated_genes_T <- DEG_Trem2 %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(300)

downregulated_genes_T <- DEG_Trem2 %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  head(200)

upregulated_genes_TNK <- DEG_TNK %>%
  filter(avg_log2FC > 0, p_val_adj <0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(100)

downregulated_genes_TNK <- DEG_TNK %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  head(100)

# Prepare gene list for enrichment
OmMet_up <- upregulated_genes
OmMet_up$gene <- rownames(upregulated_genes)

AscMet_up <- downregulated_genes
AscMet_up$gene <- rownames(downregulated_genes)

Trem2_up <- upregulated_genes_T
Trem2_up$gene <- rownames(upregulated_genes_T)

TNK_up <- upregulated_genes_TNK
TNK_up$gene <- rownames(upregulated_genes_TNK)

Asc_TNK_up <- downregulated_genes_TNK
Asc_TNK_up$gene <- rownames(downregulated_genes_TNK)

Asc_mac_up <- upregulated_genes_AM
Asc_mac_up$gene <- rownames(upregulated_genes_AM)

# Validate symbols and convert to Entrez IDs
valid_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")

filtered_OmMet <- intersect(OmMet_up$gene, valid_symbols)
entrez_ids_OmMet <- bitr(filtered_OmMet, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_OmMet <- na.omit(entrez_ids_OmMet$ENTREZID)

filtered_AscMet <- intersect(AscMet_up$gene, valid_symbols)
entrez_ids_AscMet <- bitr(filtered_AscMet, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_AscMet <- na.omit(entrez_ids_AscMet$ENTREZID)

filtered_Asc_mac_up <- intersect(Asc_mac_up$gene, valid_symbols)
entrez_ids_Asc_mac_up <- bitr(filtered_Asc_mac_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_Asc_mac_up <- na.omit(entrez_ids_Asc_mac_up$ENTREZID)

filtered_Trem2 <- intersect(Trem2_up$gene, valid_symbols)
entrez_ids_Trem2 <- bitr(filtered_Trem2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_Trem2 <- na.omit(entrez_ids_Trem2$ENTREZID)

filtered_TNK <- intersect(TNK_up$gene, valid_symbols)
entrez_ids_TNK <- bitr(filtered_TNK, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_TNK <- na.omit(entrez_ids_TNK$ENTREZID)

filtered_Asc_TNK <- intersect(Asc_TNK_up$gene, valid_symbols)
entrez_ids_Asc_TNK <- bitr(filtered_Asc_TNK, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_Asc_TNK <- na.omit(entrez_ids_Asc_TNK$ENTREZID)

### GO Biological Process Enrichment
OmMet_up_go_bp <- enrichGO(
  gene = entrez_OmMet,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(OmMet_up_go_bp, showCategory = 10, title = "Top GO Biological Processes Pathways (OmMet)")

AscMet_up_go_bp <- enrichGO(
  gene = entrez_AscMet,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(AscMet_up_go_bp, showCategory = 10, title = "Top GO Biological Processes Pathways (AscMet)")

Asc_mac_up_go_bp <- enrichGO(
  gene = entrez_Asc_mac_up,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Asc_mac_up_go_bp, showCategory = 10, title = "Top GO Biological Processes Pathways (Asc_mac)")

Trem2_up_go_bp <- enrichGO(
  gene = entrez_Trem2,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Trem2_up_go_bp, showCategory = 10, title = "Top GO Biological Processes Pathways (Trem2)")

TNK_up_go_bp <- enrichGO(
  gene = entrez_TNK,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(TNK_up_go_bp, showCategory = 10, title = "Top GO Biological Processes Pathways (TNK)")

### GO Molecular Function Enrichment
OmMet_up_go_mf <- enrichGO(
  gene = entrez_OmMet,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(OmMet_up_go_mf, showCategory = 10, title = "Top GO Molecular Function Pathways (OmMet)")

AscMet_up_go_mf <- enrichGO(
  gene = entrez_AscMet,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(AscMet_up_go_mf, showCategory = 10, title = "Top GO Molecular Function Pathways (AscMet)")

Asc_mac_up_go_mf <- enrichGO(
  gene = entrez_Asc_mac_up,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Asc_mac_up_go_mf, showCategory = 10, title = "Top GO Molecular Function Pathways (Asc MAC)")

Trem2_up_go_mf <- enrichGO(
  gene = entrez_Trem2,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Trem2_up_go_mf, showCategory = 10, title = "Top GO Molecular Function Pathways (Trem2)")

TNK_up_go_mf <- enrichGO(
  gene = entrez_TNK,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(TNK_up_go_mf, showCategory = 10, title = "Top GO Molecular Function Pathways (TNK)")

### KEGG Pathway Enrichment
OmMet_up_kegg <- enrichKEGG(
  gene = entrez_OmMet,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(OmMet_up_kegg, showCategory = 10, title = "Top KEGG Pathways (OmMet)")

AscMet_up_kegg <- enrichKEGG(
  gene = entrez_AscMet,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(AscMet_up_kegg, showCategory = 10, title = "Top KEGG Pathways (AscMet)")

Asc_mac_up_kegg <- enrichKEGG(
  gene = entrez_Asc_mac_up,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Asc_mac_up_kegg, showCategory = 10, title = "Top KEGG Pathways (Asc MAC)")

Trem2_up_kegg <- enrichKEGG(
  gene = entrez_Trem2,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Trem2_up_kegg, showCategory = 10, title = "Top KEGG Pathways (Trem2)")

TNK_up_kegg <- enrichKEGG(
  gene = entrez_TNK,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(TNK_up_kegg, showCategory = 10, title = "Top KEGG Pathways (TNK)")

### Reactome Pathway Enrichment
OmMet_up_reactome <- enrichPathway(
  gene = entrez_OmMet,
  organism = "mouse",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(OmMet_up_reactome, showCategory = 10, title = "Top Reactome Pathways (OmMet)")

AscMet_up_reactome <- enrichPathway(
  gene = entrez_AscMet,
  organism = "mouse",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(AscMet_up_reactome, showCategory = 10, title = "Top Reactome Pathways (AscMet)")

Asc_mac_up_reactome <- enrichPathway(
  gene = entrez_Asc_mac_up,
  organism = "mouse",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Asc_mac_up_reactome, showCategory = 10, title = "Top Reactome Pathways (Asc MAC)")

Trem2_up_reactome <- enrichPathway(
  gene = entrez_Trem2,
  organism = "mouse",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(Trem2_up_reactome, showCategory = 10, title = "Top Reactome Pathways (Trem2)")

TNK_up_reactome <- enrichPathway(
  gene = entrez_TNK,
  organism = "mouse",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dotplot(TNK_up_reactome, showCategory = 10, title = "Top Reactome Pathways (TNK)")

### HALLMARK ###
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# Get Hallmark gene sets for mouse (or human, adjust species if needed)
# 1. Load Hallmark gene sets for mouse
hallmark_df <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# 2. Convert your gene list from ENTREZ to SYMBOL (as you're doing)
AscMet_up_genes <- bitr(entrez_AscMet, fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Mm.eg.db)$SYMBOL

Asc_mac_up_genes <- bitr(entrez_Asc_mac_up, fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Mm.eg.db)$SYMBOL

# 3. Run enrichment analysis with SYMBOLs
AscMet_up_hallmark <- enricher(
  gene = AscMet_up_genes,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

Asc_mac_up_hallmark <- enricher(
  gene = Asc_mac_up_genes,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 4. Plot the top pathways
dotplot(Asc_mac_up_hallmark, showCategory = 10, title = "Top Hallmark Pathways (Asc MAC)")

OmMet_up_genes <- bitr(entrez_OmMet, fromType = "ENTREZID", 
                        toType = "SYMBOL", OrgDb = org.Mm.eg.db)$SYMBOL
OmMet_up_hallmark <- enricher(
  gene = OmMet_up_genes,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 4. Plot the top pathways
dotplot(OmMet_up_hallmark, showCategory = 10, title = "Top Hallmark Pathways (OmMet)")
