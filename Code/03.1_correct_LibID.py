import pandas as pd
from itertools import zip_longest
from tqdm import tqdm
import argparse
import logging

def hamming_distance(s1, s2):
    """Compute the Hamming distance between two strings, increasing the threshold for 'N' characters."""
    n_count = s1.count('N') + s2.count('N')
    base_distance = sum(c1 != c2 for c1, c2 in zip_longest(s1, s2) if c1 != 'N' and c2 != 'N')
    return base_distance, n_count

def correct_lib_id(lib_id, whitelist, base_threshold=2):
    """Correct Lib_ID based on the whitelist if Hamming distance (adjusted for 'N') is within the threshold."""
    for wl_id in whitelist:
        base_distance, n_count = hamming_distance(lib_id, wl_id)
        adjusted_threshold = base_threshold + n_count
        if base_distance < adjusted_threshold:
            return wl_id
    return lib_id

def process_files(input_file, whitelist_file, output_file, chunk_size=1000000):
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    logging.info("Loading whitelist of Lib_IDs")
    # Read the whitelist of Lib_IDs from the file
    with open(whitelist_file) as f:
        whitelist = f.read().splitlines()

    # Initialize progress bar and total row count
    logging.info("Counting total number of rows in the input file")
    total_rows = sum(1 for _ in open(input_file))
    logging.info(f"Total number of rows: {total_rows}")
    
    progress_bar = tqdm(total=total_rows, desc="Processing", unit="rows")

    # Process the input file in chunks and write to output incrementally
    logging.info("Starting to process the file in chunks")
    with pd.read_csv(input_file, sep=' ', header=None, chunksize=chunk_size, engine='python') as reader:
        for chunk in reader:
            if 2 not in chunk.columns:
                logging.error("Column index 2 does not exist in the chunk")
                break

            logging.info(f"Processing a chunk of {len(chunk)} rows")
            chunk[2] = chunk[2].apply(lambda x: correct_lib_id(x, whitelist))
            chunk.to_csv(output_file, sep=' ', header=False, index=False, mode='a')
            progress_bar.update(len(chunk))

    progress_bar.close()
    logging.info(f"Corrected data has been saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct potential sequencing errors in Lib_ID column based on a whitelist.")
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file name')
    parser.add_argument('-w', '--whitelist', type=str, required=True, help='Whitelist (Lib_ID) file name')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file name')

    args = parser.parse_args()

    process_files(args.input, args.whitelist, args.output)
