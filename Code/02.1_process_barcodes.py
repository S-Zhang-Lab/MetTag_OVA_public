import pandas as pd
import numpy as np
from itertools import combinations
import argparse
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering

# Function to calculate Hamming distance
def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

# Function to collapse barcodes based on clustering
def collapse_barcodes(barcodes, threshold):
    # Calculate pairwise Hamming distances
    distances = np.array([[hamming_distance(seq1, seq2) for seq2 in barcodes] for seq1 in barcodes])
    
    # Perform agglomerative clustering
    clustering = AgglomerativeClustering(n_clusters=None, linkage='single', distance_threshold=threshold, metric='precomputed')
    labels = clustering.fit_predict(distances)
    
    # Collapse sequences within each cluster
    collapsed_barcodes = []
    for label in np.unique(labels):
        cluster_indices = np.where(labels == label)[0]
        if len(cluster_indices) > 1:
            # Choose a representative sequence (first one in the cluster)
            representative = barcodes[cluster_indices[0]]
            collapsed_barcodes.append(representative)
        else:
            collapsed_barcodes.append(barcodes[cluster_indices[0]])
    
    return collapsed_barcodes

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Select, collapse, and merge barcode sequences based on Hamming distance threshold.')
parser.add_argument('input_file', type=str, help='Input file containing barcode sequences')
parser.add_argument('threshold', type=int, help='Hamming distance threshold for selecting and collapsing barcodes')
parser.add_argument('output_subtracted', type=str, help='Output file for barcodes not within the threshold')
parser.add_argument('output_subseted', type=str, help='Output file for barcodes within the threshold')
parser.add_argument('output_final', type=str, help='Output file for the final merged barcodes')
args = parser.parse_args()

# Read the file
data = pd.read_csv(args.input_file, sep=' ', header=None, names=['frequency', 'barcode'], dtype={'barcode': str})

# Get all barcode sequences, drop any NaNs
barcodes = data['barcode'].dropna().values

# Identify barcodes within the Hamming distance threshold
selected_indices = set()
for i, seq1 in enumerate(barcodes):
    for j, seq2 in enumerate(barcodes):
        if i < j:
            h_dist = hamming_distance(seq1, seq2)
            if h_dist < args.threshold:
                selected_indices.add(i)
                selected_indices.add(j)

# Split into selected and remaining barcodes
selected_barcodes = [barcodes[i] for i in selected_indices]
remaining_barcodes = [barcodes[i] for i in range(len(barcodes)) if i not in selected_indices]

# Save the remaining barcodes to the specified output file
remaining_df = data.iloc[[i for i in range(len(data)) if i not in selected_indices]]
remaining_df.to_csv(args.output_subtracted, index=False, header=False, sep=' ')

# Save the selected barcodes to the specified output file
selected_df = data.iloc[list(selected_indices)]
selected_df.to_csv(args.output_subseted, index=False, header=False, sep=' ')

# Collapse selected barcodes based on the Hamming distance threshold
collapsed_barcodes = collapse_barcodes(selected_barcodes, args.threshold)

# Merge the collapsed barcodes with the remaining barcodes
merged_barcodes = remaining_barcodes + collapsed_barcodes

# Save the final merged barcodes to the specified output file
final_df = pd.DataFrame(merged_barcodes, columns=['barcode'])
final_df.to_csv(args.output_final, index=False, header=False)

print(f"Remaining barcodes saved to '{args.output_subtracted}'")
print(f"Selected barcodes saved to '{args.output_subseted}'")
print(f"Final merged barcodes saved to '{args.output_final}'")
