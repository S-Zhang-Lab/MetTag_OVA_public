import argparse
import numpy as np
from annoy import AnnoyIndex
import logging
import time
import multiprocessing
import os

def setup_logging(log_filename):
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )

def hamming_distance(str1, str2):
    """Calculate the Hamming distance between two strings."""
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def barcode_to_vector(barcode):
    """Convert a barcode string to a binary vector. Assumes barcodes consist of 'A', 'C', 'G', 'T'."""
    encoding = {'A': [1, 0, 0, 0, 0], 'C': [0, 1, 0, 0, 0], 'G': [0, 0, 1, 0, 0], 'T': [0, 0, 0, 1, 0], 'N': [0, 0, 0, 0, 1]}
    return [bit for char in barcode for bit in encoding[char]]

def load_whitelist(filename):
    """Load barcodes from a whitelist file."""
    logging.info(f"Loading whitelist from {filename}")
    with open(filename) as f:
        whitelist = [line.strip() for line in f]
    logging.info(f"Loaded {len(whitelist)} barcodes from whitelist")
    return whitelist

def build_annoy_index(barcodes, vector_length, num_trees):
    """Build an Annoy index for the given barcodes."""
    logging.info(f"Building Annoy index with {num_trees} trees")
    index = AnnoyIndex(vector_length, 'hamming')
    for i, barcode in enumerate(barcodes):
        vector = barcode_to_vector(barcode)
        index.add_item(i, vector)
    
    index.build(num_trees)  # Number of trees specified by user
    logging.info("Annoy index built")
    return index

def correct_barcodes(input_list, whitelist, vector_length, hamming_threshold, column_index, num_trees, id_x):
    """Correct barcodes in the input file against the whitelist using Annoy index."""
    logging.info(f"Correcting barcodes for chunk {id_x}")
    whitelist_index = build_annoy_index(whitelist, vector_length, num_trees)
    total_barcodes = 0
    corrected_barcodes = 0
    uncorrected_barcodes = []
    final_results = []

    for ii in range(len(input_list)):
        fields = input_list[ii]
        barcode = fields[column_index]
        total_barcodes += 1
        if total_barcodes % 10000 == 0:
            logging.info(f"Processed {total_barcodes} barcodes in chunk {id_x}")

        vector = barcode_to_vector(barcode)
        nearest_indices = whitelist_index.get_nns_by_vector(vector, 50)  # Increased nearest neighbors
        corrected = False
        for idx in nearest_indices:
            nearest_barcode = whitelist[idx]
            n_number = barcode.count('N')
            if hamming_distance(barcode, nearest_barcode) - n_number <= hamming_threshold:
                fields[column_index] = nearest_barcode
                corrected = True
                corrected_barcodes += 1
                break
        if not corrected:
            fields[column_index] = barcode
            uncorrected_barcodes.append([ii, barcode])

        final_results.append(fields)

    logging.info(f"Chunk {id_x} correction completed: {corrected_barcodes} corrected, {len(uncorrected_barcodes)} uncorrected")
    return [id_x, final_results, uncorrected_barcodes]

def divide_chunks(l, n): 
    for i in range(0, len(l), n):  
        yield l[i:i + n]

def main(input_file, whitelist_file, column_index, max_distance, max_cores, num_trees):
    logging.info(f"Starting barcode correction for {input_file}")
    file = open(input_file, 'r')
    lines = file.readlines()
    for ii in range(len(lines)):
        line = lines[ii].strip().split()
        lines[ii] = line

    cpu_number = multiprocessing.cpu_count()
    if max_cores > cpu_number:
        max_cores = cpu_number
    logging.info(f"Using {max_cores} cores to process the data")

    interval = int(np.ceil(len(lines) / max_cores))
    infile = list(divide_chunks(lines, interval))
    p = multiprocessing.Pool(processes=max_cores)
    whitelist = load_whitelist(whitelist_file)

    vector_length = len(barcode_to_vector(whitelist[0]))

    results = [p.apply_async(correct_barcodes, args=(x, whitelist, vector_length, max_distance, column_index, num_trees, id_x)) for id_x, x in enumerate(infile)]
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Collapse barcodes and generate sparse matrix.")
    parser.add_argument("--input_file", default='smaller_LARRY_filtered_2percent.txt', help="Path to the input file.")
    parser.add_argument("--column_index", type=int, default=0, help="The column index for the barcodes.")
    parser.add_argument("--corrected_file", default='c_cell.txt', help="Path to the output corrected file.")
    parser.add_argument("--uncorrected_file", default='uc_cell.txt', help="Path to the output uncorrected file.")
    parser.add_argument("--max_distance", type=int, default=2, help="Maximum distance for collapsing barcodes.")
    parser.add_argument("--whitelist_file", default='CellBC_WL.txt', help="Path to the whitelist file for collapsed barcodes.")
    parser.add_argument("--max_cores", type=int, default=100, help="Maximum number of cores for parallel computation.")
    parser.add_argument("--num_trees", type=int, default=50, help="Number of trees to use in Annoy index.")
    args = parser.parse_args()

    log_filename = f"{os.path.splitext(os.path.basename(args.input_file))[0]}_barcode_correction.log"
    setup_logging(log_filename)

    t1 = time.time()
    results = main(args.input_file, args.whitelist_file, args.column_index, args.max_distance, args.max_cores, args.num_trees)
    [result.wait() for result in results]
    t2 = time.time()
    logging.info(f"Processing completed in {t2 - t1} seconds")
    logging.info('Writing results to files...')

    cor_file = open(args.corrected_file, 'w')
    uncor_file = open(args.uncorrected_file, 'w')
    pre_length = 0
    for ri in range(len(results)):
        result = results[ri].get()[1]
        for ci in range(len(result)):
            cor_file.write(' '.join(result[ci]) + '\n')

        result = results[ri].get()[2]
        for ci in range(len(result)):
            index = pre_length + result[ci][0]
            barcode = result[ci][1]
            uncor_file.write(f"{index} {barcode}\n")
        pre_length += len(results[ri].get()[1])

    cor_file.close()
    uncor_file.close()
    logging.info("Results written to files")
