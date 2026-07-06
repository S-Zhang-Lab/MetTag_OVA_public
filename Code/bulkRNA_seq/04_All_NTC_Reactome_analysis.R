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
library(stringr)
library(forcats)

setwd("/Users/s438978/Desktop/RData/RNA-seq")

# =========================================
# 1️⃣ Define count files
# =========================================
count_files <- list(
  ULA_NTC_PBS_1   = "ULA_NP_1_counts.txt",
  ULA_NTC_PBS_2   = "ULA_NP_2_counts.txt",
  ULA_NTC_PBS_3   = "ULA_NP_3_counts.txt",
  ULA_NTC_IFNg_1  = "ULA_NI_1_counts.txt",
  ULA_NTC_IFNg_2  = "ULA_NI_2_counts.txt",
  ULA_NTC_IFNg_3  = "ULA_NI_3_counts.txt",
  Reg_NTC_IFNg_1  = "Reg_NI_1_counts.txt",
  Reg_NTC_IFNg_2  = "Reg_NI_2_counts.txt",
  Reg_NTC_IFNg_3  = "Reg_NI_3_counts.txt"
)

# =========================================
# 2️⃣ Read counts
# =========================================
count_list <- lapply(count_files, function(f) {
  dt <- fread(f, skip = 1)
  as.numeric(dt[[ncol(dt)]])
})

counts <- do.call(cbind, count_list)
colnames(counts) <- names(count_files)
gene_ids <- fread(count_files[[1]], skip = 1)[[1]]
rownames(counts) <- gene_ids
counts <- apply(counts, 2, as.numeric)
rownames(counts) <- gene_ids

# -----------------------------------------
# Filter lowly expressed genes (>10 counts in at least 1 sample)
keep <- rowSums(counts > 10) >= 1
counts <- counts[keep, ]

# =========================================
# 3️⃣ Sample metadata
# =========================================
sample_info <- data.frame(
  sample = colnames(counts),
  group  = c(rep("ULA_PBS", 3),
             rep("ULA_IFNg", 3),
             rep("REG_IFNg", 3))
)
rownames(sample_info) <- sample_info$sample

# =========================================
# 4️⃣ Create DESeq2 object and normalize
# =========================================
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = sample_info,
                              design    = ~ group)
dds <- DESeq(dds)

# =========================================
# 5️⃣ Define comparisons
# =========================================
comparisons <- list(
  "ULA_IFNg_vs_PBS"      = c("group", "ULA_IFNg", "ULA_PBS"),
  "ULA_IFNg_vs_REG_IFNg" = c("group", "ULA_IFNg", "REG_IFNg")
)

# =========================================
# 6️⃣ Run DESeq2, shrink LFC, export DEGs
# =========================================
deg_list <- list()
lfc_list <- list()

for (comp_name in names(comparisons)) {
  contrast_vec <- comparisons[[comp_name]]
  
  res <- results(dds, contrast = contrast_vec)
  res <- lfcShrink(dds, contrast = contrast_vec, res = res, type = "normal")
  
  # Assign rownames if missing
  if(length(rownames(res)) == 0) rownames(res) <- rownames(counts)
  
  # Filter significant DEGs
  sig <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
  deg_list[[comp_name]] <- sig
  lfc_list[[comp_name]] <- res
  
  # Export DEGs
  if(length(sig) > 0) {
    write.csv(as.data.frame(res[sig,]), paste0("DEGs_", comp_name, ".csv"), row.names = TRUE)
  }
}

# =========================================
# 7️⃣ Full Reactome workflow (function)
# =========================================
run_reactome <- function(resLFC, direction_label="UP", padj_cutoff=0.05, log2fc_cutoff=0, top_n=15) {
  # Clean ENSEMBL IDs
  resLFC$ENSEMBL_CLEAN <- sub("\\..*$","",rownames(resLFC))
  resLFC$ENSEMBL_CLEAN <- trimws(resLFC$ENSEMBL_CLEAN)
  
  # Identify up/down genes
  if(direction_label=="UP") {
    genes <- rownames(resLFC)[resLFC$padj < padj_cutoff & resLFC$log2FoldChange > log2fc_cutoff]
  } else {
    genes <- rownames(resLFC)[resLFC$padj < padj_cutoff & resLFC$log2FoldChange < -log2fc_cutoff]
  }
  if(length(genes)==0) return(NULL)
  
  # Map ENSEMBL -> ENTREZ + SYMBOL
  valid_ids <- intersect(resLFC$ENSEMBL_CLEAN, keys(org.Mm.eg.db, keytype="ENSEMBL"))
  map_ens <- bitr(valid_ids, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Mm.eg.db)
  map_ens$log2FC <- resLFC$log2FoldChange[match(map_ens$ENSEMBL, resLFC$ENSEMBL_CLEAN)]
  lfc_by_entrez <- map_ens %>%
    filter(!is.na(ENTREZID)) %>%
    group_by(ENTREZID) %>%
    summarize(
      SYMBOLs = paste(unique(na.omit(SYMBOL)), collapse=";"),
      log2FC_gene = mean(log2FC, na.rm=TRUE),
      .groups="drop"
    )
  
  # Convert input genes to ENTREZ for Reactome
  ensembl2entrez <- bitr(genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  entrez_genes <- ensembl2entrez$ENTREZID
  if(length(entrez_genes)==0) return(NULL)
  
  # Reactome enrichment
  reactome_res <- enrichPathway(gene=entrez_genes, organism="mouse",
                                pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05,
                                readable=TRUE)
  reactome_df <- as.data.frame(reactome_res)
  if(nrow(reactome_df)==0) return(NULL)
  
  # Compute mean log2FC per pathway
  symbol_to_entrez <- bitr(unique(unlist(strsplit(as.character(reactome_df$geneID),"/"))),
                           fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  reactome_df2 <- reactome_df %>%
    rowwise() %>%
    mutate(
      genes_raw = list(strsplit(as.character(geneID),"/")[[1]]),
      genes_entrez = list(symbol_to_entrez$ENTREZID[match(genes_raw, symbol_to_entrez$SYMBOL)]),
      genes_entrez = list(genes_entrez[[1]][!is.na(genes_entrez[[1]])]),
      mean_log2FC = {
        present <- lfc_by_entrez$log2FC_gene[lfc_by_entrez$ENTREZID %in% genes_entrez[[1]]]
        if(length(present)==0) NA_real_ else mean(present, na.rm=TRUE)
      }
    ) %>%
    ungroup()
  
  # Top pathways
  plot_df <- reactome_df2 %>%
    filter(!is.na(mean_log2FC)) %>%
    arrange(if(direction_label=="UP") desc(mean_log2FC) else mean_log2FC) %>%
    slice_head(n=top_n) %>%
    mutate(Description=factor(Description, levels=rev(Description)))
  
  # Plot
  p <- ggplot(plot_df, aes(x=Description, y=mean_log2FC, fill=-log10(p.adjust))) +
    geom_col(width=0.7) +
    coord_flip() +
    scale_fill_gradient(name="-log10(adj.p)", low=if(direction_label=="UP") "#ffcccc" else "gray90",
                        high=if(direction_label=="UP") "#990000" else "gray10") +
    labs(title=paste("Reactome Pathways", direction_label, "regulated"),
         x="Pathway", y="Mean log2 Fold Change") +
    theme_minimal(base_size=13) +
    theme(axis.text.y=element_text(size=10), legend.position="right")
  
  return(list(df=reactome_df2, plot=p))
}

# =========================================
# 8️⃣ Run Reactome for all comparisons
# =========================================
reactome_results <- list()
for(comp_name in names(lfc_list)) {
  resLFC <- lfc_list[[comp_name]]
  
  # Assign rownames if missing
  if(length(rownames(resLFC))==0) rownames(resLFC) <- rownames(counts)
  
  reactome_up <- run_reactome(resLFC, direction_label="UP")
  reactome_down <- run_reactome(resLFC, direction_label="DOWN")
  
  reactome_results[[comp_name]] <- list(up=reactome_up, down=reactome_down)
  
  # Optional: display plots
  if(!is.null(reactome_up)) print(reactome_up$plot)
  if(!is.null(reactome_down)) print(reactome_down$plot)
}

# Loop over comparisons and save plots
for(comp_name in names(reactome_results)) {
  res_list <- reactome_results[[comp_name]]
  
  if(!is.null(res_list$up)) {
    pdf(paste0("Reactome_", comp_name, "_UP.pdf"), width=10, height=6)
    print(res_list$up$plot)
    dev.off()
  }
  
  if(!is.null(res_list$down)) {
    pdf(paste0("Reactome_", comp_name, "_DOWN.pdf"), width=10, height=6)
    print(res_list$down$plot)
    dev.off()
  }
}

#Methods – Gene Expression Visualization and Statistical Analysis
#Normalized gene expression values for St6gal1 were obtained from DESeq2 size-factor–adjusted counts. For visualization, expression values were plotted using box-and-whisker plots, where the box represents the interquartile range (IQR), the central line indicates the median, whiskers extend to 1.5 × IQR, and individual points correspond to biological replicates (n = 3 per condition).
#Statistical comparisons between two conditions (e.g., Reg_NTC_PBS vs ULA_NTC_PBS) were performed using a two-sample Welch’s t-test applied to normalized log-scaled expression values. When multiple pairwise comparisons were performed, p-values were adjusted using the Benjamini–Hochberg method.
#Expression values in the figure are displayed on a log10 scale to improve visual comparability across conditions.
