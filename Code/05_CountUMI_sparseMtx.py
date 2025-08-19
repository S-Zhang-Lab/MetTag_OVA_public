import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import argparse
import os

def main(input_file, output_file, csv_threshold, bc_id_file):
    # Load the data into a DataFrame
    print("Loading data...")
    data = pd.read_csv(input_file, sep=' ', header=None, names=['CellBarcode', 'UMI', 'Lib_ID', 'LARRYBarcode'])
    print("Data loaded successfully.")
    
    # Load BC_ID mapping
    print("Loading BC_ID mapping...")
    bc_id_mapping = pd.read_csv(bc_id_file, sep='\t', header=None, names=['BC_ID', 'Sequence'])
    bc_id_mapping_dict = dict(zip(bc_id_mapping['BC_ID'], bc_id_mapping['Sequence']))
    print("BC_ID mapping loaded successfully.")

    # Collapse LARRY barcodes by counting unique UMIs
    print("Collapsing LARRY barcode data...")
    collapsed_data = data.groupby(['CellBarcode', 'LARRYBarcode'])['UMI'].nunique().reset_index()
    collapsed_data.rename(columns={'UMI': 'UMICount'}, inplace=True)
    print(f"Collapsed LARRY barcode data shape: {collapsed_data.shape}")

    # Create a mapping for unique LARRY barcodes to short names
    unique_larry_barcodes = collapsed_data['LARRYBarcode'].unique()
    larry_barcode_map = {barcode: f'BC_{i+1}' for i, barcode in enumerate(unique_larry_barcodes)}
    print(f"Number of unique LARRY barcodes: {len(unique_larry_barcodes)}")

    # Map LARRY barcodes to their short names in the collapsed data
    collapsed_data['LARRYBarcodeShort'] = collapsed_data['LARRYBarcode'].map(larry_barcode_map)

    # Create a sparse matrix for LARRY barcodes
    cell_barcodes = collapsed_data['CellBarcode'].unique()
    short_larry_barcodes = collapsed_data['LARRYBarcodeShort'].unique()
    print(f"Number of unique cell barcodes: {len(cell_barcodes)}")
    print(f"Number of short LARRY barcodes: {len(short_larry_barcodes)}")

    # Create mappings from barcodes to matrix indices
    cell_barcode_to_idx = {barcode: idx for idx, barcode in enumerate(cell_barcodes)}
    short_larry_barcode_to_idx = {barcode: idx for idx, barcode in enumerate(short_larry_barcodes)}

    # Initialize the sparse matrix
    row_indices = collapsed_data['CellBarcode'].map(cell_barcode_to_idx).values
    col_indices = collapsed_data['LARRYBarcodeShort'].map(short_larry_barcode_to_idx).values
    counts = collapsed_data['UMICount'].values
    print(f"Row indices length: {len(row_indices)}")
    print(f"Col indices length: {len(col_indices)}")
    print(f"Counts length: {len(counts)}")

    # Create the sparse matrix
    matrix = sp.csr_matrix((counts, (row_indices, col_indices)), shape=(len(cell_barcodes), len(short_larry_barcodes)))

    # Save the matrix in Matrix Market format
    output_mtx_file = os.path.splitext(output_file)[0] + '.mtx'
    sio.mmwrite(output_mtx_file, matrix)
    print(f"Matrix saved as Matrix Market file: {output_mtx_file}")

    # Save cell barcodes and LARRY barcodes
    cell_barcodes_file = os.path.splitext(output_file)[0] + '_barcodes.tsv'
    with open(cell_barcodes_file, 'w') as f:
        for barcode in cell_barcodes:
            f.write(f"{barcode}\n")
    print(f"Cell barcodes saved as: {cell_barcodes_file}")

    larry_barcodes_file = os.path.splitext(output_file)[0] + '_features.tsv'
    with open(larry_barcodes_file, 'w') as f:
        for barcode in short_larry_barcodes:
            f.write(f"{barcode}\n")
    print(f"LARRY barcodes saved as: {larry_barcodes_file}")

    # Save the LARRY barcode mapping to a CSV file
    larry_mapping_df = pd.DataFrame(list(larry_barcode_map.items()), columns=['LARRYBarcode', 'ShortName'])
    mapping_csv_file = os.path.splitext(output_file)[0] + '_larry_mapping.csv'
    larry_mapping_df.to_csv(mapping_csv_file, index=False)
    print(f"LARRY barcode mapping saved as: {mapping_csv_file}")

    # Handle Lib_ID collapsing
    print("Collapsing Lib_ID data...")
    collapsed_lib_data = data.groupby(['CellBarcode', 'Lib_ID'])['UMI'].nunique().reset_index()
    collapsed_lib_data.rename(columns={'UMI': 'UMICount'}, inplace=True)
    print(f"Collapsed Lib_ID data shape: {collapsed_lib_data.shape}")

    # Create a mapping for unique Lib_IDs to short names
    unique_lib_ids = collapsed_lib_data['Lib_ID'].unique()
    lib_id_map = {lib_id: f'LIB_{i+1}' for i, lib_id in enumerate(unique_lib_ids)}
    print(f"Number of unique Lib_IDs: {len(unique_lib_ids)}")

    # Map Lib_IDs to their short names in the collapsed data
    collapsed_lib_data['LibIDShort'] = collapsed_lib_data['Lib_ID'].map(lib_id_map)

    # Create a sparse matrix for Lib_IDs
    short_lib_ids = collapsed_lib_data['LibIDShort'].unique()
    print(f"Number of short Lib_IDs: {len(short_lib_ids)}")

    # Create mappings from Lib_IDs to matrix indices
    short_lib_id_to_idx = {lib_id: idx for idx, lib_id in enumerate(short_lib_ids)}

    # Initialize the sparse matrix
    lib_row_indices = collapsed_lib_data['CellBarcode'].map(cell_barcode_to_idx).values
    lib_col_indices = collapsed_lib_data['LibIDShort'].map(short_lib_id_to_idx).values
    lib_counts = collapsed_lib_data['UMICount'].values
    print(f"Lib row indices length: {len(lib_row_indices)}")
    print(f"Lib col indices length: {len(lib_col_indices)}")
    print(f"Lib counts length: {len(lib_counts)}")

    # Create the sparse matrix for Lib_IDs
    lib_matrix = sp.csr_matrix((lib_counts, (lib_row_indices, lib_col_indices)), shape=(len(cell_barcodes), len(short_lib_ids)))

    # Save the Lib_ID matrix in Matrix Market format
    output_lib_mtx_file = os.path.splitext(output_file)[0] + '_lib.mtx'
    sio.mmwrite(output_lib_mtx_file, lib_matrix)
    print(f"Lib_ID matrix saved as Matrix Market file: {output_lib_mtx_file}")

    # Save the Lib_ID barcode mapping to a CSV file
    lib_mapping_df = pd.DataFrame(list(lib_id_map.items()), columns=['Lib_ID', 'ShortName'])
    lib_mapping_csv_file = os.path.splitext(output_file)[0] + '_lib_mapping.csv'
    lib_mapping_df.to_csv(lib_mapping_csv_file, index=False)
    print(f"Lib_ID mapping saved as: {lib_mapping_csv_file}")

    # Save the Lib_ID barcodes as a TSV file
    lib_id_barcodes_file = os.path.splitext(output_file)[0] + '_lib_id_barcode.tsv'
    with open(lib_id_barcodes_file, 'w') as f:
        for barcode in unique_lib_ids:
            f.write(f"{barcode}\n")
    print(f"Lib_ID barcodes saved as: {lib_id_barcodes_file}")

    # Save to CSV if the size is small
    if len(collapsed_data) <= csv_threshold:
        csv_output_file = os.path.splitext(output_file)[0] + '.csv'
        # Use short LARRY barcodes in the CSV output
        collapsed_data_short = collapsed_data[['CellBarcode', 'LARRYBarcodeShort', 'UMICount']].copy()
        collapsed_data_short.rename(columns={'LARRYBarcodeShort': 'LARRYBarcode'}, inplace=True)
        collapsed_data_short.to_csv(csv_output_file, index=False)
        print(f"Data is small, saved as CSV: {csv_output_file}")
    else:
        print(f"Data is large, not saved as CSV.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process DNA barcode data and save as sparse matrix.')
    parser.add_argument('input_file', type=str, help='Path to the input file containing the DNA barcode data.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    parser.add_argument('bc_id_file', type=str, help='Path to the BC_ID mapping file.')
    parser.add_argument('--csv_threshold', type=int, default=100000, help='Row count threshold for saving as CSV (default: 100000).')

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.csv_threshold, args.bc_id_file)
