### Step7.3.2 analysis UTSW02 - cancer cell focused ###

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

#Load in transcriptome R object and check metadata columns

Ova_merged_SCT <- readRDS("Ova_merged_SCT_EA_061024.rds")
head(Ova_merged_SCT@meta.data)
DimPlot(Ova_merged_SCT, label = TRUE)
### after 07.3.2 preprocessing it will be: Ova_merged_SCT <- readRDS("Ova_merged_SCT_EA_061824.rds")
### can also load: Mets <- readRDS("Mets_SCT_EA_062024.rds") directly

### 0 - start with updating the metadata with LARRY BC assignments from SZ object
#Load in SZs data and examine structure
Ova_merged_SCT_SZ <- readRDS("MetTag.ova_L1L2_merged_Post_SCT_integrated.rds")
head(Ova_merged_SCT_SZ@meta.data)
updated.meta <- Ova_merged_SCT_SZ@meta.data
table(updated.meta$Combined.LARRY_group)
table(updated.meta$om_LibID_first_bc)
table(updated.meta$Asc_LibID_first_bc)

# Extract the columns
asc_first_bc <- updated.meta$Asc_LARRY_first_bc
om_first_bc <- updated.meta$om_LARRY_first_bc
# Find common values
common_values <- intersect(asc_first_bc, om_first_bc)
# Display the common values
print(common_values)

# Combine existing and new metadata
# Extract original clusters with cell barcodes -- if I don't do this first, clustering will change
original_clusters <- data.frame(
  cell_barcodes = rownames(Ova_merged_SCT@meta.data),
  seurat_clusters = Ova_merged_SCT@meta.data$seurat_clusters
)

# Merge the original clusters with updated metadata based on cell barcodes
updated.meta.2 <- merge(
  updated.meta,
  original_clusters,
  by.x = "row.names",
  by.y = "cell_barcodes",
  all.x = TRUE
)

# Set the row names to the cell barcodes
rownames(updated.meta.2) <- updated.meta.2$Row.names
updated.meta.2$Row.names <- NULL

Ova_merged_SCT <- AddMetaData(object = Ova_merged_SCT, metadata = updated.meta.2)

# Check if clustering remained the same after integrating metadata 
Idents(Ova_merged_SCT) <- "seurat_clusters"
DimPlot(Ova_merged_SCT, label = TRUE)

# Now subset cancer cell clusters only into a new Mets object
Idents(Ova_merged_SCT) <- "Combined.HTO_group"
Met.samples <- subset(Ova_merged_SCT, idents = c("asc:1_AscMet", "om:2_OmMet"))
Idents(Met.samples) <- "seurat_clusters"
Mets <- subset(Met.samples, idents = c("0", "1", "4", "10", "11", "24", "26", "38"))

# Set seed
set.seed(42)
#re-cluster subsetted cell type
Mets <- RunPCA(Mets, verbose = FALSE)
Idents(Mets) <- "seurat_clusters"
Mets <- RunUMAP(Mets, dims = 1:10, verbose = FALSE)
Mets <- FindNeighbors(Mets, dims = 1:30, verbose = FALSE)
Mets <- FindClusters(Mets, verbose = FALSE, resolution = 0.2)
DimPlot(Mets, label = TRUE)
DimPlot(Mets, group.by = "Combined.HTO_group")
#export dimplots

### 01 - Examine differences between ascites and omentum cancer cell RNA
#Compare Cluster Proportions Between Control and Experimental 
table(Mets@meta.data$seurat_clusters,Mets@meta.data$Combined.HTO_group)
freq_table <- prop.table(x = table(Mets@meta.data$seurat_clusters,Mets@meta.data$Combined.HTO_group),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(Mets@meta.data$seurat_clusters) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)

#Examine Differentially Expressed Genes Among RNA Clusters
Idents(Mets) <- "seurat_clusters"
Mets <- PrepSCTFindMarkers(Mets)
Mets.markers <- FindAllMarkers(Mets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Mets.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
View(Mets.markers)
top10_Mets_DEG <- Mets.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Mets, features = top10_Mets_DEG$gene) +scale_fill_gradientn(colors = c("blue", "white", "red")) 
write.csv(Mets.markers,"Mets_marker_genes_subsetted.csv")

VlnPlot(Mets, features = c("Cxcl1", "Mki67", "Col6a1", "Cd274", "B2m", "Hbb-bs"))

# Cancer cell subtypes (not very precise classification)
# Cluster 0 - Chemokine high
# Cluster 1 - Collagen high
# Cluster 2 - Proliferative
# Cluster 3 - Hemoglobin gene high
# Cluster 4 - PD-L1 high MHC-1 low
# Cluster 5 - Hemoglobin gene high
# Cluster 6 - PD-L1 low MHC-1 high

# Examine DEGs between Omental metastases and ascites
Idents(Mets) <- "Combined.HTO_group"
DEG_Mets <- FindMarkers(Mets, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_Mets, file = "OmMet_compared_to_ascMet_DEG.csv")

DEG_Mets$neg_adj.pVal <- -1*log(DEG_Mets$p_val_adj, 10)
pThresh <- 10
DEG_Mets$hits <- ifelse(DEG_Mets$hits > pThresh, c("hits"), c("")) 
DEG_Mets$label <- ifelse(DEG_Mets$hits == "hits", row.names(DEG_Mets), c("")) 
plot(DEG_Mets$avg_log2FC, DEG_Mets$neg_adj.pVal, col=ifelse(DEG_Mets$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 300), pch=ifelse(DEG_Mets$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-10, 10, by=2.5), col="lightgray")
abline(h=seq(0,300, by=50), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Mets$avg_log2FC, DEG_Mets$neg_adj.pVal, col=ifelse(DEG_Mets$neg_adj.pVal > pThresh, ifelse(DEG_Mets$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Mets$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Mets$neg_adj.pVal > pThresh, 1.4, 1.2))
# export volcano 6x5

ggplotly(ggplot(data = DEG_Mets,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Mets), colour=ifelse(DEG_Mets$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

### 02_Perform similar DEG analysis but with BC.ID filtering
# Compare highly enriched IDs vs low

# Check transcriptome differences among BC.IDs within individual treatment groups

Idents(Mets) <- "om_LibID_first_bc"
DEG_Om.Mets.BCID6.vs.BCID12 <- FindMarkers(Mets, ident.1 = c("R1_BC_ID6"), ident.2 = c("R1_BC_ID1", "R1_BC_ID2"), verbose = TRUE, logfc.threshold = 0.05)
Idents(Mets) <- "Asc_LibID_first_bc"
DEG_Asc.Mets.BCID6.vs.BCID12 <- FindMarkers(Mets, ident.1 = c("R1_BC_ID6"), ident.2 = c("R1_BC_ID1", "R1_BC_ID2"), verbose = TRUE, logfc.threshold = 0.05)
# no DEGs for Asc comparison, but there are DEGs in Om - add to CRISPR list
write.csv(DEG_Om.Mets.BCID6.vs.BCID12, file = "DEG_Om.Mets.BCID6.vs.BCID12.csv")

### Qunantify BC.ID frequency across each HTO in ascmet and ommet for statistical analysis in Prism Graphpad (converting from counts to frequency in Excel first)
# Load required packages
library(dplyr)
library(tidyr)

# Subset the data
Idents(Mets) <- "Combined.HTO_group"
OmMet <- subset(Mets, idents = c("om:2_OmMet"))
AscMet <- subset(Mets, idents = c("asc:1_AscMet"))

# Extract metadata
OmMet.meta <- OmMet@meta.data
AscMet.meta <- AscMet@meta.data

# Function to create frequency table for each Combined.HTO_group
create_frequency_table <- function(meta_data, group_column, barcode_column) {
  meta_data %>%
    group_by(across(all_of(group_column)), across(all_of(barcode_column))) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = all_of(barcode_column), values_from = count, values_fill = list(count = 0))
}

# Create frequency tables
OmMet_freq_table <- create_frequency_table(OmMet.meta, "om.HTO_classification", "om_LibID_first_bc")
AscMet_freq_table <- create_frequency_table(AscMet.meta, "Asc.HTO_classification", "Asc_LibID_first_bc")

# Print the frequency tables
print(OmMet_freq_table)
print(AscMet_freq_table)

# make volcano plot for Om BC.ID6 DEGs:

DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal <- -1*log(DEG_Om.Mets.BCID6.vs.BCID12$p_val_adj, 10)
pThresh <- 2
plot(DEG_Om.Mets.BCID6.vs.BCID12$avg_log2FC, DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal, col=ifelse(DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 7), pch=ifelse(DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-7.5, 7.5, by=2.5), col="lightgray")
abline(h=seq(0,7, by=1.5), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Om.Mets.BCID6.vs.BCID12$avg_log2FC, DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal, col=ifelse(DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal > pThresh, ifelse(DEG_Om.Mets.BCID6.vs.BCID12$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal > pThresh, 1.4, 1.2))
# export volcano 6x5

ggplotly(ggplot(data = DEG_Om.Mets.BCID6.vs.BCID12,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Om.Mets.BCID6.vs.BCID12), colour=ifelse(DEG_Om.Mets.BCID6.vs.BCID12$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))


### 03_ Examine clonality information - Top BC compared to others

Idents(Mets) <- "Combined.HTO_group"
OmMet <- subset(Mets, idents = c("om:2_OmMet"))
OmMet@meta.data$BC_status <- ifelse(OmMet@meta.data$Combined.LARRY_group != "not_detected-om", "LARRY_detected", "LARRY_not_detected")
Idents(OmMet) <- "BC_status"
OmMet <- subset(OmMet, idents = c("LARRY_detected"))
Idents(OmMet) <- "Combined.LARRY_group"

# Draw a frequency map using ggplot
Om_id_freq <- table(OmMet$Combined.LARRY_group)
Om_freq_df <- data.frame(MULTI_ID = names(Om_id_freq), Frequency = as.integer(Om_id_freq))
Om_freq_df_sorted <- Om_freq_df[order(-Om_freq_df$Frequency), ]
Om_freq_df_sorted$MULTI_ID <- factor(Om_freq_df_sorted$MULTI_ID, levels = Om_freq_df_sorted$MULTI_ID[order(-Om_freq_df_sorted$Frequency)])
Om_freq_df_sorted <- Om_freq_df_sorted[-1, ]
#above code removes the first row which contains "LARRY_ratio_less 2 -om" frequencies instead of actual BCs

ggplot(Om_freq_df_sorted, aes(x = MULTI_ID, y = Frequency)) + geom_bar(stat = "identity", fill = "blue") + theme_minimal() + labs(title = "Frequency of BCs", x = "MULTI_ID", y = "Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot only top 61 (more aesthetically pleasing) -- matches all unique clones in AscMet 
top_61_df <- head(Om_freq_df_sorted, 61)
ggplot(top_61_df, aes(x = MULTI_ID, y = Frequency)) + geom_bar(stat = "identity", fill = "blue") + theme_minimal() + labs(title = "Frequency of BCs", x = "MULTI_ID", y = "Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

Idents(Mets) <- "Combined.HTO_group"
AscMet <- subset(Mets, idents = c("asc:1_AscMet"))
AscMet@meta.data$BC_status <- ifelse(AscMet@meta.data$Combined.LARRY_group != "not_detected-asc", "LARRY_detected", "LARRY_not_detected")
Idents(AscMet) <- "BC_status"
AscMet <- subset(AscMet, idents = c("LARRY_detected"))
Idents(AscMet) <- "Combined.LARRY_group"

# Draw a frequency map using ggplot
Asc_id_freq <- table(AscMet$Combined.LARRY_group)
Asc_freq_df <- data.frame(MULTI_ID = names(Asc_id_freq), Frequency = as.integer(Asc_id_freq))
Asc_freq_df_sorted <- Asc_freq_df[order(-Asc_freq_df$Frequency), ]
Asc_freq_df_sorted$MULTI_ID <- factor(Asc_freq_df_sorted$MULTI_ID, levels = Asc_freq_df_sorted$MULTI_ID[order(-Asc_freq_df_sorted$Frequency)])
Asc_freq_df_sorted <- Asc_freq_df_sorted[-1, ]
#above code removes the first row which contains "LARRY_ratio_less 2 -asc" frequencies instead of actual BCs

ggplot(Asc_freq_df_sorted, aes(x = MULTI_ID, y = Frequency)) + geom_bar(stat = "identity", fill = "orange") + theme_minimal() + labs(title = "Frequency of BCs", x = "MULTI_ID", y = "Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Take a closer look at top BCs compared to others in OmMets
# First: create a new metadata slot called Expansion value
Om.meta <- OmMet@meta.data
barcode_counts <- table(Om.meta$Combined.LARRY_group)
Om.meta$Expansion_value <- "Mid"
expanded_barcodes <- names(barcode_counts[barcode_counts >= 10])
non_expanded_barcodes <- names(barcode_counts[barcode_counts == 1])
not_detected_barcodes <- "LARRY_ratio_less 2-om"

Om.meta$Expansion_value[Om.meta$Combined.LARRY_group %in% expanded_barcodes] <- "Expanded"
Om.meta$Expansion_value[Om.meta$Combined.LARRY_group %in% non_expanded_barcodes] <- "Non-expanded"
Om.meta$Expansion_value[Om.meta$Combined.LARRY_group %in% not_detected_barcodes] <- "Not-detected"

# Add the updated metadata back to the Seurat object
OmMet <- AddMetaData(object = OmMet, metadata = Om.meta)
# Verify the new metadata slot
table(OmMet@meta.data$Expansion_value)
Idents(OmMet) <- "Expansion_value"

# Repeat same process for AscMet
Asc.meta <- AscMet@meta.data
barcode_counts <- table(Asc.meta$Combined.LARRY_group)
Asc.meta$Expansion_value <- "Mid"
expanded_barcodes <- names(barcode_counts[barcode_counts >= 3])
non_expanded_barcodes <- names(barcode_counts[barcode_counts == 1])
not_detected_barcodes <- "LARRY_ratio_less 2-asc"

Asc.meta$Expansion_value[Asc.meta$Combined.LARRY_group %in% expanded_barcodes] <- "Expanded"
Asc.meta$Expansion_value[Asc.meta$Combined.LARRY_group %in% non_expanded_barcodes] <- "Non-expanded"
Asc.meta$Expansion_value[Asc.meta$Combined.LARRY_group %in% not_detected_barcodes] <- "Not-detected"

# Add the updated metadata back to the Seurat object
AscMet <- AddMetaData(object = AscMet, metadata = Asc.meta)
# Verify the new metadata slot
table(AscMet@meta.data$Expansion_value)
Idents(AscMet) <- "Expansion_value"

## Go back to OmMets and look for DEGs between expanded and non-expanded clones
DEG_OmMet_expansion <- FindMarkers(OmMet, ident.1 = "Expanded", ident.2 = "Non-expanded", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_OmMet_expansion, file = "DEG_OmMet_expansion.csv")

# Filter for significant upregulated genes in both datasets
sig_up_OmMet_expansion <- DEG_OmMet_expansion[DEG_OmMet_expansion$p_val_adj < 0.05 & DEG_OmMet_expansion$avg_log2FC > 0, ]
sig_up_OmMets_BCID6_vs_BCID12 <- DEG_Om.Mets.BCID6.vs.BCID12[DEG_Om.Mets.BCID6.vs.BCID12$p_val_adj < 0.05 & DEG_Om.Mets.BCID6.vs.BCID12$avg_log2FC > 0, ]

# Load the libraries
library(DESeq2)
library(VennDiagram)

# Find the common upregulated genes
common_genes <- intersect(rownames(sig_up_OmMet_expansion), rownames(sig_up_OmMets_BCID6_vs_BCID12))
# Extract the common genes from the first dataset (you can choose either)
common_genes_data <- sig_up_OmMet_expansion[common_genes, ]
# Write to CSV
write.csv(common_genes_data, "common_upregulated_genes.csv")

# Plot the Venn diagram
venn.plot <- venn.diagram(
  x = list(OmMet_expansion = rownames(sig_up_OmMet_expansion), OmMets_BCID6_vs_BCID12 = rownames(sig_up_OmMets_BCID6_vs_BCID12)),
  category.names = c("OmMet_expansion", "OmMets_BCID6_vs_BCID12"),
  filename = NULL,
  output = TRUE
)

# Save the Venn diagram
png("Venn_diagram.png")
grid.draw(venn.plot)
dev.off()

# make volano plot for clonality comparison - later will do GSVA
DEG_OmMet_expansion$neg_adj.pVal <- -1*log(DEG_OmMet_expansion$p_val_adj, 10)
pThresh <- 2
plot(DEG_OmMet_expansion$avg_log2FC, DEG_OmMet_expansion$neg_adj.pVal, col=ifelse(DEG_OmMet_expansion$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 10), pch=ifelse(DEG_OmMet_expansion$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-6, 6, by=2), col="lightgray")
abline(h=seq(0,10, by=2), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_OmMet_expansion$avg_log2FC, DEG_OmMet_expansion$neg_adj.pVal, col=ifelse(DEG_OmMet_expansion$neg_adj.pVal > pThresh, ifelse(DEG_OmMet_expansion$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_OmMet_expansion$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_OmMet_expansion$neg_adj.pVal > pThresh, 1.4, 1.2))

#make the web version to add gene names later
ggplotly(ggplot(data = DEG_OmMet_expansion,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_OmMet_expansion), colour=ifelse(DEG_OmMet_expansion$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

# Verify that the top clone is mainly enriched in one HTO (I did not save this graph)
Idents(OmMet) <- "Combined.LARRY_group"
OmMet.BC4 <- subset(OmMet, idents = c("BC-4-om"))
OmMet.BC49 <- subset(OmMet, idents = c("BC-49-om"))
OmMet.BC6 <- subset(OmMet, idents = c("BC-6-om"))
p1 <- DimPlot(OmMet.BC4)
p2 <- DimPlot(OmMet.BC49)
p3 <- DimPlot(OmMet.BC6)


### 04 - Generate heatmap showing enrichment of clones across HTOs (repeat same for Asc samples)
library(tidyr)
library(pheatmap)
Om.meta <- OmMet@meta.data
Om.meta <- Om.meta %>% filter(Combined.LARRY_group != "LARRY_ratio_less 2-om")
barcode_summary <- Om.meta %>% group_by(om.HTO_classification, Combined.LARRY_group) %>% summarize(count = n()) %>% pivot_wider(names_from = om.HTO_classification, values_from = count, values_fill = list(count = 0))

Asc.meta <- AscMet@meta.data
Asc.meta <- Asc.meta %>% filter(Combined.LARRY_group != "LARRY_ratio_less 2-asc")
barcode_summary <- Asc.meta %>% group_by(Asc.HTO_classification, Combined.LARRY_group) %>% summarize(count = n()) %>% pivot_wider(names_from = Asc.HTO_classification, values_from = count, values_fill = list(count = 0))

#Calculate the total count for each barcode across all HTO samples
barcode_summary$total_count <- rowSums(barcode_summary[, -1])

# Ensure barcode_summary is a data.frame
barcode_summary <- as.data.frame(barcode_summary)

# Remove row names if they exist
rownames(barcode_summary) <- NULL

# Convert columns to appropriate types
barcode_summary$Combined.LARRY_group <- as.character(barcode_summary$Combined.LARRY_group)
barcode_summary$total_count <- as.numeric(barcode_summary$total_count)

# Step 1: Order the data frame by total_count in descending order
ordered_summary <- barcode_summary[order(-barcode_summary$total_count), ]

# Step 2: Subset the top 61 rows
top_61_summary <- ordered_summary[1:61, ]

# Step 3: Extract the Combined.LARRY_group column
top_barcodes <- top_61_summary$Combined.LARRY_group

# Check the result
print(top_barcodes)

# Remove the total_count column
top_61_summary <- top_61_summary %>%
  select(-total_count)

barcode_summary_normalized <- top_61_summary %>% mutate(across(-Combined.LARRY_group, ~ ./ sum(.)))

barcode_matrix <- as.matrix(barcode_summary_normalized[, -1])
rownames(barcode_matrix) <- barcode_summary_normalized$Combined.LARRY_group
custom_colors <- colorRampPalette(c("white", "blue", "purple"))(50)
# Create the heatmap
pheatmap(barcode_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE,
         color = custom_colors,
         main = "LARRY Barcodes Enrichment Across HTO Samples")

saveRDS(Mets, file = "Mets_SCT_EA_062024.rds")
saveRDS(Ova_merged_SCT, file = "Ova_Merged_SCT_EA_061824.rds")


### 05_Examine DEGs BETWEEN top3 expanded clones
Idents(OmMet) <- "Combined.LARRY_group"

DEG_BC4_BC49 <- FindMarkers(OmMet, ident.1 = "BC-4-om", ident.2 = "BC-49-om", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_BC4_BC49, file = "DEG_BC4_BC49.csv")

DEG_BC4_BC1116 <- FindMarkers(OmMet, ident.1 = "BC-4-om", ident.2 = "BC-1116-om", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_BC4_BC1116, file = "DEG_BC4_BC1116.csv")

DEG_BC4_BC6 <- FindMarkers(OmMet, ident.1 = "BC-4-om", ident.2 = "BC-6-om", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_BC4_BC6, file = "DEG_BC4_BC6.csv")

DEG_BC49_BC6 <- FindMarkers(OmMet, ident.1 = "BC-49-om", ident.2 = "BC-6-om", verbose = TRUE, logfc.threshold = 0.05)
write.csv(DEG_BC49_BC6, file = "DEG_BC49_BC6.csv")

# After OvCa CRISPR screen completion - check if Gbp2b is a DEG in BCID6 Asc specifically
# Subsetting cells from "om:2_OmMet" group that belong to "R1_BC_ID6"
om_R1_BC_ID6 <- subset(Mets, subset = Combined.HTO_group == "om:2_OmMet" & om_LibID_first_bc == "R1_BC_ID6")

# Subsetting cells from "asc:1_AscMet" group that belong to "R1_BC_ID6"
asc_R1_BC_ID6 <- subset(Mets, subset = Combined.HTO_group == "asc:1_AscMet" & Asc_LibID_first_bc == "R1_BC_ID6")

# Merging the two subsets
merged_R1_BC_ID6 <- merge(om_R1_BC_ID6, y = asc_R1_BC_ID6)
Idents(merged_R1_BC_ID6) <- "Combined.HTO_group"
merged_R1_BC_ID6 <- PrepSCTFindMarkers(merged_R1_BC_ID6)
# Perform the DEG analysis between "om:2_OmMet" and "asc:1_AscMet"
DEGs_R1_BC_ID6 <- FindMarkers(merged_R1_BC_ID6, ident.1 = "om:2_OmMet", ident.2 = "asc:1_AscMet", min.pct = 0.25)
write.csv(DEGs_R1_BC_ID6, file = "DEG_OmMet_vs_AscMet_BCID6.csv")

### 06_Showcase dormancy program in ascites
# Extract metadata
OmMet.meta <- subset(Mets, idents = "om:2_OmMet")@meta.data
AscMet.meta <- subset(Mets, idents = "asc:1_AscMet")@meta.data

# Count frequency of each barcode
library(dplyr)

OmMet_summary <- OmMet.meta %>%
  group_by(om_LARRY_first_bc) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  rename(Combined.LARRY_group = om_LARRY_first_bc) %>%
  mutate(Site = "Omentum")

AscMet_summary <- AscMet.meta %>%
  group_by(Asc_LARRY_first_bc) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  rename(Combined.LARRY_group = Asc_LARRY_first_bc) %>%
  mutate(Site = "Ascites")

top_n <- 50  # or any number you prefer

# Get top barcodes for each condition separately
top_om_bcs <- OmMet_summary %>%
  arrange(desc(total_count)) %>%
  slice_head(n = top_n)

top_asc_bcs <- AscMet_summary %>%
  arrange(desc(total_count)) %>%
  slice_head(n = top_n)

combined_top_bc <- bind_rows(top_om_bcs, top_asc_bcs)

combined_top_bc <- combined_top_bc %>%
  mutate(Freq_Bin = cut(
    total_count,
    breaks = c(-Inf, 1, 5, Inf),
    labels = c("Freq. 1", "Freq. 2–5", "Freq. 5+"),
    right = TRUE
  ))

# Remove barcodes labeled "not_detected"
combined_top_bc <- combined_top_bc %>%
  filter(!grepl("not_detected", Combined.LARRY_group))

library(ggplot2)

# Summarize for plotting
bc_bin_summary <- combined_top_bc %>%
  group_by(Site, Freq_Bin) %>%
  summarise(Count = n(), .groups = "drop")

# Plot proportions
ggplot(bc_bin_summary, aes(x = Site, y = Count, fill = Freq_Bin)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(
    values = c(
      "Freq. 1"   = "gray40",
      "Freq. 2–5" = "gray80",
      "Freq. 5+"  = "red"
    )
  ) +
  ylab("Proportion of Barcodes") +
  ggtitle("Barcode Expansion Categories (Top 50)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



### 07
# Load necessary libraries for plotting
library(EnhancedVolcano)

# Generating a volcano plot from DEGs
EnhancedVolcano(DEGs_R1_BC_ID6,
                lab = rownames(DEGs_R1_BC_ID6),
                x = 'avg_log2FC',  # Log fold change
                y = 'p_val_adj',  # Adjusted p-value
                title = 'DEGs between om:2_OmMet and asc:1_AscMet (R1_BC_ID6)',
                pCutoff = 0.05,  # Significance threshold
                FCcutoff = 0.25,
                selectLab = c("Gbp2b", "Marco")
)

# Enhanced Volcano with manual gene selection and adjustments
genes_to_label <- c("Gbp2b", "Marco", "Ms4a8a", "Cd300lb", "Gm46224")  # Replace with your actual gene names

EnhancedVolcano(DEGs_R1_BC_ID6,
                lab = rownames(DEGs_R1_BC_ID6),            # Gene labels
                x = 'avg_log2FC',                              # Log2 fold change column
                y = 'p_val_adj',                           # Adjusted p-value column
                title = 'DEGs between OmMet and AscMet (BC.ID6)',
                pCutoff = 0.05,                            # p-value cutoff
                FCcutoff = 0.25,                           # Fold change cutoff
                pointSize = 3.0,                           # Size of points
                labSize = 3.5,                             # Size of gene labels
                selectLab = genes_to_label,                # Manually select genes to label
                max.overlaps = Inf,                        # Prevent label overlap suppression
                drawConnectors = TRUE,                     # Draw lines connecting labels to points
                widthConnectors = 0.5,                     # Width of the lines
                colConnectors = 'grey30',                  # Color of the connector lines
                xlim = c(-10, 10),                           # Adjust x-axis limits for visibility
                ylim = c(0, 100))  # Adjust y-axis limits for visibility

# Display top3 clones onto a UMAP
# Set seed
set.seed(42)
#re-cluster subsetted cell type OmMet
Idents(OmMet) <- "seurat_clusters"
OmMet <- RunPCA(OmMet, verbose = FALSE)
OmMet <- RunUMAP(OmMet, dims = 1:10, verbose = FALSE)
OmMet <- FindNeighbors(OmMet, dims = 1:30, verbose = FALSE)
OmMet <- FindClusters(OmMet, verbose = FALSE, resolution = 0.2)
DimPlot(OmMet, label = TRUE)

#re-cluster subsetted cell type AscMet
Idents(AscMet) <- "seurat_clusters"
AscMet <- RunPCA(AscMet, verbose = FALSE)
AscMet <- RunUMAP(AscMet, dims = 1:10, verbose = FALSE)
AscMet <- FindNeighbors(AscMet, dims = 1:30, verbose = FALSE)
AscMet <- FindClusters(AscMet, verbose = FALSE, resolution = 0.2)
DimPlot(AscMet, label = TRUE)
# export OmMet and AscMet RDS


# Subset the Seurat object to include only the cells in the 3 specified groups
subset_cells <- subset(OmMet, Combined.LARRY_group %in% c("BC-4-om", "BC-1116-om", "BC-6-om"))
# Set the identities of the cells to the 'Combined.LARRY_group'
Idents(subset_cells) <- "Combined.LARRY_group"
# Plot the UMAP with different colors for the three cell populations
DimPlot(subset_cells, reduction = "umap", group.by = "Combined.LARRY_group", 
        cols = c("BC-4-om" = "red", "BC-1116-om" = "blue", "BC-6-om" = "green"),
        pt.size = 0.5) +
  labs(title = "BC-4-om, BC-1116-om, and BC-6-om cells")

# Create a frequency table for Combined.LARRY_group and seurat_clusters
freq_table <- table(subset_cells$Combined.LARRY_group, subset_cells$seurat_clusters)

# Convert the table to a data frame for easy plotting
freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("Population", "SeuratCluster", "Frequency")

# Create a frequency table for Combined.LARRY_group and seurat_clusters
freq_table <- table(subset_cells$Combined.LARRY_group, subset_cells$seurat_clusters)

# Load ggplot2
library(ggplot2)

# Plot the bar chart
ggplot(freq_df, aes(x = Population, y = Frequency, fill = SeuratCluster)) +
  geom_bar(stat = "identity", position = "dodge") +  # Create a grouped bar chart
  labs(title = "Frequency of Seurat Clusters by Population", 
       x = "Population", 
       y = "Frequency") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Optional: Adjust color palette

Idents(subset_cells) <- "seurat_clusters"
# Find DEGs between cluster 0 and all other clusters
DEGs_cluster_0_vs_others <- FindMarkers(subset_cells, ident.1 = "0", ident.2 = c("1", "2", "3", "4"))

# Find all markers OmMet & AscMet seperate RDS objects
OmMet.markers <- FindAllMarkers(OmMet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top30 <- OmMet.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DoHeatmap(OmMet, features = top30$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

AscMet.markers <- FindAllMarkers(AscMet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top30 <- AscMet.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DoHeatmap(AscMet, features = top30$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))
# save heatmaps as PDF

# Check Gbp2b levels (CRISPR screen results) + proliferation markers
FeaturePlot(Mets, features = c("Gbp2b", "Marco", "Cd300lb", "Ms4a8a", "Slfn1"), reduction = "umap")

# Violin plot comparing expression across ascites and omentum groups
# Step 1: Set identity to Combined.HTO_group
Idents(Mets) <- "Combined.HTO_group"

# Step 2: Subset each group
AscMet <- subset(Mets, idents = "asc:1_AscMet")
OmMet <- subset(Mets, idents = "om:2_OmMet")

# Step 3: Downsample OmMet to match AscMet cell number
set.seed(42)  # for reproducibility
n_asc <- ncol(AscMet)
OmMet_downsampled <- subset(OmMet, cells = sample(colnames(OmMet), n_asc))

# Step 4: Merge downsampled datasets
combined <- merge(AscMet, OmMet_downsampled)

# Step 5: Violin plots with p-values using Seurat + ggpubr
library(Seurat)
library(ggplot2)
library(ggpubr)

# Genes of interest
genes <- c("Mki67", "Top2a", "Ifngr1")

# Step 6: Generate violin plots and p-values
plots <- lapply(genes, function(gene) {
  VlnPlot(combined, features = gene, group.by = "Combined.HTO_group", pt.size = 0.1) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
})

# Step 7: Arrange and display plots
library(patchwork)
wrap_plots(plots, ncol = 1)

# Create a data frame to hold p-values
pval_df <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Loop through genes, calculate p-values
for (gene in genes) {
  expr <- FetchData(combined, vars = c(gene, "Combined.HTO_group"))
  colnames(expr) <- c("expression", "group")
  
  # Wilcoxon rank-sum test
  test <- wilcox.test(expression ~ group, data = expr)
  
  # Save to data frame
  pval_df <- rbind(pval_df, data.frame(Gene = gene, P_Value = test$p.value))
}

# Print the p-values
print(pval_df)
#Gene    P_Value
#1  Mki67 0.06555212
#2  Top2a 0.04242041
#3 Ifngr1 0.03450140

# Plot function with red dot for mean
plot_with_mean <- function(gene) {
  p <- VlnPlot(combined, features = gene, group.by = "Combined.HTO_group", pt.size = 0.1) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "red", color = "black") +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

# Generate plots
plots <- lapply(genes, plot_with_mean)

# Display vertically
library(patchwork)
wrap_plots(plots, ncol = 1)


### 06_GATING out Gbp2b high cells ###
# reload Mets object if necessary
Mets <- readRDS("Mets_SCT_EA_062024.rds")
# Make two sub-populations: Gbp2b high and Gbp2b low
Idents(Mets) <- "seurat_clusters"
plotStat1xGbp2b <- FeatureScatter(Mets, feature1 = "Stat1", feature2 = "Gbp2b")
#select Stat1+/- only, no Gbp2b expression
Mets <- CellSelector(plot = plotStat1xGbp2b, object = Mets, ident = "Gbp2b(low)")
#select Gbp2b high cells
Mets <- CellSelector(plot = plotStat1xGbp2b, object = Mets, ident = "Gbp2b(high)")
DimPlot(Mets, split.by = "Combined.HTO_group")
# export 4 x 4 pdf of gate
DEG_Gbp2bhigh <- FindMarkers(Mets, ident.1 = "Gbp2b(high)", ident.2 = "Gbp2b(low)", verbose = TRUE, logfc.threshold = 0.05)
write.csv(as.matrix(DEG_Gbp2bhigh), file = "DEG_Gbp2bHigh_vs_Low.csv")

# make volcano plot - later will do GSVA
DEG_Gbp2bhigh$neg_adj.pVal <- -1*log(DEG_Gbp2bhigh$p_val_adj, 10)
pThresh <- 10
plot(DEG_Gbp2bhigh$avg_log2FC, DEG_Gbp2bhigh$neg_adj.pVal, col=ifelse(DEG_Gbp2bhigh$neg_adj.pVal > pThresh,"red3","black"), ylim=c(0, 300), pch=ifelse(DEG_Gbp2bhigh$neg_adj.pVal > pThresh, 19, 1), cex=1.2)
abline(v=seq(-5, 15, by=2.5), col="lightgray")
abline(h=seq(0,300, by=50), col="lightgray")
abline(h=pThresh, col="black", lty=2)
points(DEG_Gbp2bhigh$avg_log2FC, DEG_Gbp2bhigh$neg_adj.pVal, col=ifelse(DEG_Gbp2bhigh$neg_adj.pVal > pThresh, ifelse(DEG_Gbp2bhigh$avg_log2FC>0,"red3","dodgerblue3"),"black"), pch=ifelse(DEG_Gbp2bhigh$neg_adj.pVal > pThresh, 19, 1), cex=ifelse(DEG_Gbp2bhigh$neg_adj.pVal > pThresh, 1.4, 1.2))

#make the web version to add gene names later
ggplotly(ggplot(data = DEG_Gbp2bhigh,aes(x=avg_log2FC, y=neg_adj.pVal ,text=rownames(DEG_Gbp2bhigh), colour=ifelse(DEG_Gbp2bhigh$neg_adj.pVal > pThresh,"darkblue","darkred")),cex=0.5) + geom_point(shape=21, size=3, fill="white") + geom_hline(yintercept = 0) + theme(panel.background = element_rect(fill=NA),panel.grid.major = element_line(colour = "grey80"),panel.ontop = TRUE))

# Extract cell identities (Gbp2b High vs Low)
Mets$Gbp2b_status <- Idents(Mets)

# Check how many cells in each group
table(Mets$Gbp2b_status)

# Extract metadata from Mets (subsetted cancer cell object)
gbp2b_metadata <- Mets@meta.data[, "Gbp2b_status", drop=FALSE]

# Add it to full dataset (cells not in Mets will be NA)
Ova_merged_SCT <- readRDS("Ova_merged_SCT_EA_061824.rds")
Ova_merged_SCT <- AddMetaData(Ova_merged_SCT, metadata = gbp2b_metadata)

### 07 Check CRISPR hits via violin plot + FOLR1 expression on cancer cells
Mets@meta.data$Combined.HTO_group.order <- recode(Mets@meta.data$Combined.HTO_group, "asc:1_AscMet" = "01_AscMet", "om:2_OmMet" = "02_OmMet")
Idents(Mets) <- "Combined.HTO_group"
VlnPlot(Mets, features = c("Ccr1", "Cd274"))

# Specify the order of the groups for plotting
group_order <- c("01_AscMet", "02_OmMet")

# Plot the VlnPlot with the specified group order and gene expression for Ifng
VlnPlot(Mets, features = c("Gbp2b", "Marco", "Ms4a8a", "Cd300lb"), group.by = "Combined.HTO_group", 
        idents = group_order, 
        pt.size = 1)

# Subset the Seurat object to only include AscMet and OmMet
Mets_subset <- subset(Mets, Combined.HTO_group %in% c("asc:1_AscMet", "om:2_OmMet"))

# Get the number of cells in the smaller group (AscMet)
n_cells_AscMet <- sum(Mets_subset$Combined.HTO_group == "asc:1_AscMet")

# Get the cells from the OmMet group and downsample to match AscMet
OmMet_cells <- WhichCells(Mets_subset, ident = "om:2_OmMet")

# Randomly sample the same number of cells as AscMet
set.seed(42)  # Set a seed for reproducibility
OmMet_cells_downsampled <- sample(OmMet_cells, size = n_cells_AscMet)

# Subset the Seurat object to only include the downsampled OmMet cells and all AscMet cells
Mets_balanced <- subset(Mets_subset, cells = c(WhichCells(Mets_subset, ident = "asc:1_AscMet"), OmMet_cells_downsampled))

# Now, plot the VlnPlot for Ifng expression
VlnPlot(Mets_balanced, features =  c("Gbp2b", "Marco", "Ms4a8a", "Cd300lb", "Slfn1"), group.by = "Combined.HTO_group.order", 
        pt.size = 1)

VlnPlot(Mets_balanced, features =  c("Ifngr1"), group.by = "Combined.HTO_group.order", pt.size = 1)

library(ggpubr)

# Create a violin plot for the gene of interest
# Check the presence of the gene in the correct assay
# Ensure you're accessing the correct assay
# Access the normalized expression data from the SCT assay
gene_expression_data <- Mets_balanced[["SCT"]]@data

# Check if the gene exists in the expression data
gene_of_interest <- "Ifngr1"  # Replace with your gene name

if (gene_of_interest %in% rownames(gene_expression_data)) {
  
  # Create the violin plot with p-value annotation
  plot <- VlnPlot(Mets_balanced, features = gene_of_interest, group.by = "Combined.HTO_group.order") +
    stat_compare_means(
      comparisons = list(c("01_AscMet", "02_OmMet")),  # Specify which groups to compare
      label = "p.value",  # Show p-value as a significance level (e.g., *, **, ***, etc.)
      size = 6,
      label.y = 5.5
    )
  
  # Print the plot
  print(plot)
  
} else {
  message("Gene not found in the dataset.")
}

# Subset only the two groups of interest
gene_expr_subset <- subset(gene_expr, group %in% c("01_AscMet", "02_OmMet"))

# Perform Wilcoxon test
p_val <- compare_means(expression ~ group, data = gene_expr_subset,
                       method = "wilcox.test")

# Print p-value to console
print(p_val)

# Results for IFNGR1
# A tibble: 1 × 8
#.y.        group1   group2        p p.adj p.format p.signif method  
#<chr>      <chr>    <chr>     <dbl> <dbl> <chr>    <chr>    <chr>   
#  1 expression 02_OmMet 01_AscMet 0.600   0.6 0.6      ns       Wilcoxon
#  
# Create the violin plot and add median points
plot <- VlnPlot(Mets_balanced, features = gene_of_interest, group.by = "Combined.HTO_group.order") +
  stat_summary(
    fun = "mean",  # Calculate the median
    geom = "point",  # Display the median as a point
    size = 3,        # Size of the point
    color = "red",    # Color of the point
    shape = 18        # Shape of the point (18 is for a solid circle)
  ) +
  stat_compare_means(
    comparisons = list(c("01_AscMet", "02_OmMet")),  # Specify which groups to compare
    label = "p.value",  # Show exact p-value
    size = 5            # Adjust size of p-value label
  )

# Print the plot
print(plot)




