# to use: python 02_Create_LARRY_WL.py --input_file filtered_LARRY.txt --output_file LARRY_WL.txt --freq_output_file LARRY_freq.txt --plot_file knee_plot.pdf --sensitivity 5
import argparse
import collections
import numpy as np
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm
from kneed import KneeLocator

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Parse command line arguments
parser = argparse.ArgumentParser(description="Process LARRY barcodes and determine knee point.")
parser.add_argument('--input_file', type=str, required=True, help='Input file name (filtered barcodes only)')
parser.add_argument('--output_file', type=str, required=True, help='Output file name (whitelist file)')
parser.add_argument('--freq_output_file', type=str, required=True, help='Barcode frequency output file name')
parser.add_argument('--plot_file', type=str, default='knee_plot.pdf', help='Knee plot PDF file name')
parser.add_argument('--sensitivity', type=float, required=True, help='Knee detection sensitivity (recommended range: 0.5 to 2.0)')
args = parser.parse_args()

# Assign command line arguments to variables
input_file = args.input_file
output_file = args.output_file
freq_output_file = args.freq_output_file
plot_file = args.plot_file
sensitivity = args.sensitivity

# Initialize a counter for the barcodes
barcode_counter = collections.Counter()

logging.info("Starting to read the data")

# Step 1: Read the data
with open(input_file, 'r') as file:
    lines = file.readlines()
    for line in tqdm(lines, desc="Processing lines"):
        barcode = line.strip()
        barcode_counter[barcode] += 1

logging.info("Completed reading the data")

# Step 2: Sort the barcodes by frequency
logging.info("Sorting the barcodes by frequency")
sorted_barcodes = sorted(barcode_counter.items(), key=lambda item: item[1], reverse=True)
frequencies = np.array([count for barcode, count in sorted_barcodes])
logging.info("Barcodes sorted")

# Step 3: Calculate the knee point using kneed
logging.info("Calculating the knee point using kneed")
x = np.arange(len(frequencies))
kneedle = KneeLocator(x, frequencies, curve='convex', direction='decreasing', S=sensitivity)
knee_point = kneedle.knee
logging.info(f"Knee point calculated at index {knee_point}")
print(f"Knee point calculated at index {knee_point} with frequency {frequencies[knee_point]}")

# Save the knee plot as a PDF file
logging.info("Generating knee plot")
plt.plot(frequencies, label='Frequency')
plt.axvline(x=knee_point, color='r', linestyle='--', label='Knee Point')
plt.xlabel('Barcode Rank')
plt.ylabel('Frequency')
plt.title('Knee Plot for LARRY Barcode Frequencies')
plt.legend()
plt.savefig(plot_file)
plt.close()
logging.info(f"Knee plot saved as {plot_file}")

# Select the barcodes up to the knee point
logging.info("Selecting the predicted real barcodes")
predicted_real_barcodes = [barcode for barcode, count in sorted_barcodes[:knee_point+1]]

# Step 4: Save the predicted real LARRY barcodes to a file
logging.info("Saving the predicted real LARRY barcodes to file")
with open(output_file, 'w') as file:
    for barcode in tqdm(predicted_real_barcodes, desc="Writing barcodes to file"):
        file.write(f"{barcode}\n")

logging.info(f"Predicted real LARRY barcodes saved to {output_file}")

# Step 5: Save the frequencies to a file
logging.info("Saving the barcode frequencies to file")
with open(freq_output_file, 'w') as file:
    for barcode, count in tqdm(sorted_barcodes, desc="Writing frequencies to file"):
        file.write(f"{count} {barcode}\n")

logging.info(f"Barcode frequencies saved to {freq_output_file}")
