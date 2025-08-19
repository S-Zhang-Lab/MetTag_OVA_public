import pandas as pd
import argparse
import logging
from tqdm import tqdm

# Set up logging
logging.basicConfig(filename='postmatching.log', level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Function to read the whitelist from a file
def read_whitelist(file_path):
    logging.info(f"Reading whitelist from {file_path}")
    with open(file_path, 'r') as file:
        whitelist = [line.strip() for line in file.readlines()]
    logging.debug(f"Whitelist contents: {whitelist[:10]}... (showing first 10 entries)")
    return whitelist

# Function to debug and clean data
def clean_and_match(series, whitelist):
    logging.info("Cleaning and matching data")
    # Remove leading/trailing spaces and make everything uppercase
    cleaned_series = series.str.strip().str.upper()
    cleaned_whitelist = [item.strip().upper() for item in whitelist]
    matches = cleaned_series.isin(cleaned_whitelist).sum()
    logging.debug(f"Matches found: {matches} out of {len(series)}")
    return matches, len(series)

def main(args):
    # Read the data file in chunks
    logging.info(f"Reading data file {args.data_file} in chunks")
    chunk_size = 100000
    total_rows = 0
    df_chunks = []

    try:
        with tqdm(total=sum(1 for _ in open(args.data_file)), desc="Loading data", unit="rows") as pbar:
            for chunk in pd.read_csv(args.data_file, header=None, delim_whitespace=True, names=['Column1', 'Column2', 'Column3', 'Column4'], chunksize=chunk_size, engine='c'):
                df_chunks.append(chunk)
                total_rows += len(chunk)
                pbar.update(len(chunk))
    except pd.errors.ParserError:
        logging.warning("C engine failed, switching to python engine")
        with tqdm(total=sum(1 for _ in open(args.data_file)), desc="Loading data", unit="rows") as pbar:
            for chunk in pd.read_csv(args.data_file, header=None, delim_whitespace=True, names=['Column1', 'Column2', 'Column3', 'Column4'], chunksize=chunk_size, engine='python'):
                df_chunks.append(chunk)
                total_rows += len(chunk)
                pbar.update(len(chunk))

    df = pd.concat(df_chunks, ignore_index=True)
    logging.info(f"Finished loading data. Total rows: {total_rows}")
    print(f"Finished loading data. Total rows: {total_rows}")

    # Read the whitelists
    logging.info("Reading whitelists")
    cellwhite_list = read_whitelist(args.cellwhite_list_file)
    larry_white_list = read_whitelist(args.larry_white_list_file)
    bc_id_white_list = read_whitelist(args.bc_id_white_list_file)

    # Count matches for column 1 (CellBarcode)
    logging.info("Counting matches for Cell_BC (Column 1)")
    print("Counting matches for Cell_BC (Column 1)")
    matches_col1, total_col1 = clean_and_match(df['Column1'], cellwhite_list)
    print(f"Matches in Cell_BC (Column 1): {matches_col1} out of {total_col1}")

    # Count matches for column 4 (LARRYBarcode)
    logging.info("Counting matches for LARRY_BC (Column 4)")
    print("Counting matches for LARRY_BC (Column 4)")
    matches_col4, total_col4 = clean_and_match(df['Column4'], larry_white_list)
    print(f"Matches in LARRY_BC (Column 4): {matches_col4} out of {total_col4}")

    # Count matches for column 3 (BC_ID)
    logging.info("Counting matches for BC_ID (Column 3)")
    print("Counting matches for BC_ID (Column 3)")
    matches_col3, total_col3 = clean_and_match(df['Column3'], bc_id_white_list)
    print(f"Matches in BC_ID (Column 3): {matches_col3} out of {total_col3}")

    # Calculate percentages
    percent_col1 = (matches_col1 / total_col1) * 100 if total_col1 > 0 else 0
    percent_col4 = (matches_col4 / total_col4) * 100 if total_col4 > 0 else 0
    percent_col3 = (matches_col3 / total_col3) * 100 if total_col3 > 0 else 0

    # Logging results
    logging.info(f"Total number of lines: {total_rows}")
    logging.info(f"Number of matches in Cell_BC: {matches_col1}")
    logging.info(f"Percentage of matches (corrected) in CellBC: {percent_col1:.2f}%")
    logging.info(f"Number of matches in LARRY_BC: {matches_col4}")
    logging.info(f"Percentage of matches (corrected) in LARRY_BC: {percent_col4:.2f}%")
    logging.info(f"Number of matches in BC_ID: {matches_col3}")
    logging.info(f"Percentage of matches (corrected) in BC_ID: {percent_col3:.2f}%")

    # Print results
    print(f"Total number of lines: {total_rows}")
    print(f"Number of matches in Cell_BC: {matches_col1}")
    print(f"Percentage of matches (corrected) in CellBC: {percent_col1:.2f}%")
    print(f"Number of matches in LARRY_BC: {matches_col4}")
    print(f"Percentage of matches (corrected) in LARRY_BC: {percent_col4:.2f}%")
    print(f"Number of matches in BC_ID: {matches_col3}")
    print(f"Percentage of matches (corrected) in BC_ID: {percent_col3:.2f}%")

    # Save results to a file
    results = {
        'Total number of lines': [total_rows],
        'Number of matches in Cell_BC': [matches_col1],
        'Percentage of matches (corrected) in CellBC': [percent_col1],
        'Number of matches in LARRY_BC': [matches_col4],
        'Percentage of matches (corrected) in LARRY_BC': [percent_col4],
        'Number of matches in BC_ID': [matches_col3],
        'Percentage of matches (corrected) in BC_ID': [percent_col3]
    }
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output_file, index=False)
    logging.info(f"Results saved to {args.output_file}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Count matches in data file against whitelists.")
    parser.add_argument('data_file', type=str, help="Path to the data file.")
    parser.add_argument('cellwhite_list_file', type=str, help="Path to the Cellwhite list (c1).")
    parser.add_argument('larry_white_list_file', type=str, help="Path to the LARRY whitelist (c4).")
    parser.add_argument('bc_id_white_list_file', type=str, help="Path to the BC_ID whitelist.")
    parser.add_argument('output_file', type=str, help="Path to the output file.")

    # Parse arguments and run main function
    args = parser.parse_args()
    main(args)
