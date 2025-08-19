### NICHENET ANALYSIS ON CANCER CELLS INTERACTING WITH MACS (R02)###
# Goal: assess differential interactions between OmMet and AscMet cancer cells
# with macs 

# load required packages
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(SeuratObject)
library(circlize)
library(RColorBrewer)
library(dplyr)



# load Seurat object of interest
setwd("C:/Users/Zhang Lab/Desktop/Emma/R01_and_2_w_LARRY")
Ova_merged_SCT <- readRDS("Ova_merged_SCT_EA_061824.rds")

seuratObj <- Ova_merged_SCT
seuratObj@meta.data %>% head()
seuratObj@meta.data$seurat_clusters %>% table()
seuratObj@meta.data$Combined.HTO_group %>% table()
seuratObj@meta.data$Broad_cell.ID %>% table()

# For older Seurat objects, you may need to run this
seuratObj <- UpdateSeuratObject(seuratObj)

# Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks
organism = "mouse"

ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final_mouse.rds")
ligand_target_matrix[1:5, 1:5] # target genes in rows, ligands in columns
lr_network <- readRDS("lr_network_mouse_21122021.rds")
head(lr_network)

weighted_networks <- readRDS("weighted_networks_nsga2r_final_mouse.rds")
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
head(weighted_networks$lr_sig)
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# Define receiver cells - cancer
# Subset the Seurat object to focus on Myeloid sender cells

### Define sender and receiver cells
Idents(seuratObj) <- "Broad_cell.ID"
receiver <- "ID8.Tumor"
sender_celltypes <- c("Myeloid")

expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)
expressed_genes_sender <- unique(unlist(lapply(sender_celltypes, get_expressed_genes, seuratObj, 0.05)))

all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)


### Define gene sets of interest to compare - AscMet vs OmMet
condition_oi <- "asc:1_AscMet"
condition_reference <- "om:2_OmMet"

seurat_obj_receiver <- subset(seuratObj, idents = receiver)
seurat_obj_receiver <- PrepSCTFindMarkers(seurat_obj_receiver)

DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,
                                 ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "Combined.HTO_group",
                                 min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

### Ligand activity prediction
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

### Visualize ligand activities
p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) + 
  geom_histogram(color = "black", fill = "darkorange") + 
  geom_vline(aes(xintercept = min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color = "red", linetype = "dashed", size = 1) + 
  labs(x = "Ligand activity (AUPR)", y = "# ligands") +
  theme_classic()

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

make_heatmap_ggplot(vis_ligand_aupr,
                    "Prioritized ligands", "Ligand activity", 
                    legend_title = "AUPR", color = "red") + 
  theme(axis.text.x.top = element_blank())


### Active target gene inference
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

top_n_genes <- active_ligand_target_links_df %>%
  arrange(desc(weight)) %>%
  head(200)

filtered_active_ligand_target_links_df <- top_n_genes

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = filtered_active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "blue", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke", high = "blue")


### Receptor plot
lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from, to)
best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df <- lr_network_top_df_large %>% spread("from", "weight", fill = 0)
lr_network_top_matrix <- lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors <- dist(lr_network_top_matrix, method = "binary")
hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
order_receptors <- hclust_receptors$labels[hclust_receptors$order]

dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()

p_ligand_receptor_network <- vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands", "Receptors", color = "darkred", x_axis_position = "top", legend_title = "Prior interaction potential")
p_ligand_receptor_network



ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  ggplot(aes(x = reorder(test_ligand, aupr_corrected), y = aupr_corrected)) +
  geom_col(fill = "darkred") +
  coord_flip() +
  labs(x = "Ligand", y = "AUPR Score", title = "Top ligands ranked by activity") +
  theme_minimal()

top_lr <- weighted_networks_lr %>%
  filter(from %in% best_upstream_ligands, to %in% expressed_receptors) %>%
  group_by(from) %>%
  top_n(3, wt = weight)

ggplot(top_lr, aes(x = from, y = to, size = weight)) +
  geom_point(color = "steelblue") +
  theme_minimal() +
  labs(x = "Ligand", y = "Receptor", size = "Interaction strength") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Condition specific ligand usage
VlnPlot(seuratObj, features = best_upstream_ligands, group.by = "Combined.HTO_group", split.by = "Broad_cell.ID")

#DotPlot for ligand expression 
DotPlot(seuratObj, features = best_upstream_ligands, group.by = "Broad_cell.ID") +
  RotatedAxis() +
  ggtitle("Ligand expression in macrophages")

library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse annotation package
library(dplyr)

# Assuming active_ligand_target_links_df is your dataframe with ligand-target links
# and best_upstream_ligands is your vector of top mouse ligands

# Filter for links from top ligands only
top_ligand_targets <- active_ligand_target_links_df %>%
  filter(ligand %in% best_upstream_ligands) %>%
  pull(target) %>%       # Extract the target genes
  unique()               # Keep unique targets only

# Run GO enrichment on these predicted target genes
ego_targets <- enrichGO(gene = top_ligand_targets,
                        OrgDb = org.Mm.eg.db,
                        keyType = "SYMBOL",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

# Visualize top enriched GO terms
dotplot(ego_targets, showCategory = 10) + 
  ggtitle("GOMF enrichment of predicted target genes")

### VENN DIAGRAM OF HUMAN & MOUSE INTERSECTION

best_upstream_ligands_mouse <- best_upstream_ligands
best_upstream_receptors_mouse <- best_upstream_receptors
expressed_genes_sender_mouse <- expressed_genes_sender
robust_ligands_mouse <- intersect(best_upstream_ligands_mouse, expressed_genes_sender_mouse)

best_upstream_ligands_human <- best_upstream_ligands
best_upstream_receptors_human <- best_upstream_receptors
expressed_genes_sender_human <- expressed_genes_sender
robust_ligands_human <- intersect(best_upstream_ligands_human, expressed_genes_sender_human)

mouse_to_human <- c(
  Bst2 = "BST2",
  Ccl12 = "CCL2",
  Cd14 = "CD14",
  Cd48 = "CD48",
  Csf1 = "CSF1",
  Edil3 = "EDIL3",
  Ifitm6 = "IFITM3",
  Il10 = "IL10",
  Il1b = "IL1B",
  Il1rn = "IL1RN",
  Il27 = "IL27",
  Jam2 = "JAM2",
  Ocln = "OCLN",
  Selp = "SELP",
  Tgfb1 = "TGFB1",
  Tnf = "TNF"
)

robust_ligands_mouse <- mouse_to_human[robust_ligands_mouse]

shared_ligands <- intersect(robust_ligands_human, robust_ligands_mouse)

ligand_sets <- list(Human = robust_ligands_human, Mouse = robust_ligands_mouse)

library(VennDiagram)

venn_lig <- venn.diagram(
  x = ligand_sets,
  category.names = c("Human", "Mouse"),
  filename = NULL,
  fill = c("red", "darkred"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2
)

grid.newpage(); grid.draw(venn_lig)
write.csv(data.frame(shared_ligands), "shared_ligands_human_mouse.csv", row.names = FALSE)
