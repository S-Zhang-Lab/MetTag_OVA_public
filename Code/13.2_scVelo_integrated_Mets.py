# Code adapted from SG - 02-19-25#
# Load required packages in RStudio
library(Seurat)
library(SeuratDisk)
library(anndata)
library(SeuratObject)

# set working directory
setwd("/Users/s438978/Desktop/RData/MetTag_U6_LARRY_BC_OVA/scVelo2_SZ")

# load rds saved from Seurat
Mets <- readRDS("Mets_only_post_clusterGESA_with_expansion_06-04-25.rds")

# filter out problematic columns that might induce errors with SaveH5Seurat
Mets@meta.data[] <- lapply(Mets@meta.data, function(x) {
  if (is.factor(x)) as.character(x) else x
})
if ("X" %in% colnames(Mets@meta.data)) {
  Mets@meta.data$X <- NULL
}
colnames(Mets@meta.data)[duplicated(colnames(Mets@meta.data))]
Mets@meta.data <- Mets@meta.data[, !duplicated(colnames(Mets@meta.data))]

# Convert Assay5 manually to Assay
convert_assay5_to_assay <- function(assay5_obj) {
  counts_matrix <- GetAssayData(assay5_obj, slot = "counts")
  data_matrix <- GetAssayData(assay5_obj, slot = "data")
  
  new_assay <- CreateAssayObject(counts = counts_matrix)
  new_assay <- SetAssayData(new_assay, slot = "data", new.data = data_matrix)
  
  return(new_assay)
}

# Apply conversion
Mets@assays <- lapply(Mets@assays, function(assay) {
  if (inherits(assay, "Assay5")) {
    return(convert_assay5_to_assay(assay))
  } else {
    return(assay)
  }
})

Mets@meta.data[] <- lapply(Mets@meta.data, function(x) {
  if (is.factor(x)) return(as.character(x))
  if (!is.numeric(x) && !is.character(x)) return(as.character(x))
  return(x)
})
Mets@meta.data <- Mets@meta.data[, !duplicated(colnames(Mets@meta.data))]

SaveH5Seurat(Mets, filename = "Mets.h5Seurat", overwrite = TRUE)
Convert("Mets.h5Seurat", dest = "h5ad") # the resulting h5ad can be loaded in python using > adata_subset = sc.read_h5ad("Mets.h5ad")

# save metadata table:
Mets$barcode <- colnames(Mets)
Mets$UMAP_1 <- Mets@reductions$umap@cell.embeddings[,1]
Mets$UMAP_2 <- Mets@reductions$umap@cell.embeddings[,2]
write.csv(Mets@meta.data, file='Mets_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(Mets, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(Mets@reductions$pca@cell.embeddings, file="pca.csv", quote=F, row.names=F)

# write gene names
write.table(
  data.frame("gene"=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F)

# >> In bash >>
# module load python/3.8.x-anaconda
# conda activate ISS # activate a specific conda enironnment

# >> If in Window Office PC terminal/pwoershell, use conda env as below: 
# cd C:\Users\s208205\Downloads\R_projects\LoomSG
# conda activate python310
# python

import os
import loompy
import scanpy as sc
import scvelo as scv
import anndata as ad
import pandas as pd
import numpy as np
import cellrank as cr

# Define the base directory where loom files are stored
base_dir = "/Users/s438978/Desktop/RData/MetTag_U6_LARRY_BC_OVA/Looms_L12"

# Function to generate the full path for a loom file
def get_loom_path(filename):
    return os.path.join(base_dir, filename)

# Dictionary containing original loom file paths
original_looms = {
    "L1": get_loom_path("R02_EA_L1_LARRYreseq.loom"),
    "L2": get_loom_path("R02_EA_L2_LARRYreseq.loom"),
}

# Define the output directory for modified loom files
output_base_dir = os.path.join(base_dir, "modified_looms")
os.makedirs(output_base_dir, exist_ok=True)  # Ensure output directory exists

# Function to modify CellIDs and create a new loom file
def add_suffix_to_loom(loom_path, original_suffix, new_suffix, output_path):
    with loompy.connect(loom_path) as ds:
        # Modify CellID: Remove the original suffix (if present) and append new suffix
        new_cell_ids = [
            id.split(':')[1].replace(original_suffix, '') + new_suffix if original_suffix in id else id + new_suffix
            for id in ds.ca["CellID"]
        ]
        
        # Create a new loom file with updated CellIDs
        loompy.create(output_path, ds.layers, ds.ra, {"CellID": new_cell_ids})

# Process each loom file
for label, loom_path in original_looms.items():
    output_path = os.path.join(output_base_dir, f"modified_{label}.loom")
    
    add_suffix_to_loom(loom_path, original_suffix="x", new_suffix=f"_{label}", output_path=output_path)
    
    print(f"Created: {output_path}")

# Original suffix to remove (assuming "x" is the character to be removed)
original_suffix = "x"

# New suffixes to append for each loom file
new_suffixes = ["_L1", "_L2"]

# Modify cell IDs in each loom file
temp_looms = []
for (label, original), new_suffix in zip(original_looms.items(), new_suffixes):
    temp_output_path = os.path.join(output_base_dir, f"modified_{label}.loom")
    add_suffix_to_loom(original, original_suffix, new_suffix, temp_output_path)
    temp_looms.append(temp_output_path)

# ✅ Check the suffixes
for temp_loom in temp_looms:
    with loompy.connect(temp_loom) as ds:
        cell_ids = ds.ca["CellID"]  # Extract modified cell IDs
    print(f"First 10 cell IDs in {temp_loom}: {cell_ids[:10]}")

# ✅ Combine modified loom files into a single loom
output_file = os.path.join(base_dir, "combined.loom")
loompy.combine(temp_looms, output_file, key="Accession")
print(f"Combined loom saved at: {output_file}")

# ✅ Check combined loom file
with loompy.connect(output_file) as ds:
    cell_ids = ds.ca["CellID"]  # Extract cell IDs
print(f"First 10 cell IDs in combined loom: {cell_ids[:10]}")

# ✅ Check for specific suffix patterns (_L1 and _L2)
for suffix in ["_L1", "_L2"]:
    count = sum(id.endswith(suffix) for id in cell_ids)
    print(f"Number of cells with suffix '{suffix}': {count}")

# ✅ Load the combined loom file into an AnnData object for scVelo
adata = sc.read_loom(output_file)
print(f"Loaded AnnData object with shape: {adata.shape}")

# Now `combined.loom` is saved and ready for scVelo analysis.


# ==========  scanpy sc analysis ===========
from scipy import io
import anndata

X = io.mmread("counts.mtx")

from anndata import AnnData

X = AnnData(X=X.transpose().tocsr())

cell_meta = pd.read_csv("Mets_metadata.csv")
with open("gene_names.csv", 'r') as f: gene_names = f.read().splitlines()
>>> print(type(X))
<class 'anndata._core.anndata.AnnData'>
>>> from anndata import AnnData
adata = AnnData(X=X.X.tocsr())  # No transpose
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
sc.pl.umap(adata, color=['Combined.HTO_group'], frameon=False, save=True)
adata.write('my_data.h5ad')
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo',facecolor='white', dpi=100, frameon = False)
cr.settings.verbosity = 2
ldata = sc.read_loom('combined.loom')

# Remove "Om_" or "Asc_" prefix
adata.obs.index = adata.obs.index.str.replace(r"^(Om_|Asc_)", "", regex=True)

# Replace _L1 and _L2 with -1
ldata.obs.index = ldata.obs.index.str.replace(r"_L[12]$", "-1", regex=True)

common_cells = adata.obs.index[adata.obs.index.isin(ldata.obs.index)]
print(f"Number of matching cells: {len(common_cells)}")

barcodes = ldata.obs.index.tolist()
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata.obs.index = barcodes
adata = scv.utils.merge(adata, ldata)

# 7. Perform RNA velocity analysis:

sc.pl.umap(adata, color='Combined.HTO_group', frameon=False, legend_loc='on data', title='', save='_HTO.pdf')
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='', save='_clusters.pdf')
scv.pl.proportions(adata, groupby='celltype_full')
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
>>> scv.tl.velocity_graph(adata)
computing velocity graph (using 1/8 cores)
  0%|          | 0/5328 [00:00<?, ?cells/s]
/Users/ealeksa2/opt/anaconda3/lib/python3.9/site-packages/scvelo/core/_parallelize.py:138: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray.
  res = np.array(res) if as_array else res
    finished (0:00:39) --> added 
    'velocity_graph', sparse matrix with cosine correlations (adata.uns)
>>> scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
computing velocity embedding
    finished (0:00:01) --> added
    'velocity_umap', embedded velocity vectors (adata.obsm)
saving figure to file ./figures/scvelo_embedding.pdf
scv.pl.velocity_embedding_grid(adata, basis='umap', color='Combined_HTO.group', save='embedding_grid_HTO.pdf', title='', scale=0.25)
>>> scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save='embedding_grid.pdf', title='', scale=0.25)
saving figure to file ./figures/scvelo_embedding_grid.pdf


# Individual genes
scv.pl.velocity(adata, ['Marco',  'Gbp2b', 'Ms4a8a', 'Cd300lb'], ncols=2)
scv.pl.velocity_graph(adata, threshold=.1)

import matplotlib.pyplot as plt
# Pseudotime analysis
scv.tl.velocity_pseudotime(adata)
fig, ax = plt.subplots(figsize=(10, 10))
# Generate scatter plot
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', ax=ax, show=False)
plt.savefig("pseudotime_plot.pdf", format="pdf")

# PAGA velocity graph
import matplotlib.pyplot as plt

scv.pl.paga(adata, basis='umap', show=False)
plt.savefig("my_paga_plot.pdf", format="pdf")




