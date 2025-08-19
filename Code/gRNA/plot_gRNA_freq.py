import pandas as pd
import matplotlib.pyplot as plt

# Prompt the user for the input file path and output PDF file name
input_file_path = input("Enter the path to the input file: ")
output_pdf_path = input("Enter the name for the output PDF file: ")

# Load the data, assuming your file handling issue is resolved with delim_whitespace
data = pd.read_csv(input_file_path, delim_whitespace=True, names=['Frequency', 'Sequence'])

# Ensure the 'Frequency' column is numeric
data['Frequency'] = pd.to_numeric(data['Frequency'])

# Filter out gRNA with fewer than 5 reads
filtered_data = data[data['Frequency'] >= 10]

# Plotting the histogram
plt.figure(figsize=(10, 6))
plt.hist(filtered_data['Frequency'], bins=100, color='blue', alpha=0.7)
plt.title('Histogram of gRNA Sequence Frequencies')
plt.xlabel('Frequency')
plt.ylabel('Count of Sequences')
plt.grid(True)

# Save the plot as a PDF file
plt.savefig(output_pdf_path, format='pdf')
