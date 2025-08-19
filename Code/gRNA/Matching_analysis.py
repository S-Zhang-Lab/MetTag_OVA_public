import pandas as pd

# Function to request file path and numeric inputs
def get_input(prompt, type_=str):
    user_input = input(prompt)
    if type_ == int:
        while not user_input.isdigit():  # Validate for integer input
            print("Please enter a valid integer.")
            user_input = input(prompt)
        return int(user_input)
    return user_input

# Load the CRISPR guide data with user input
guide_path = get_input("Enter the path to the CRISPR guide data file: ")
guides = pd.read_csv(guide_path)
print("First few rows of CRISPR guide data before cleaning:")
print(guides.head())

# Remove 'AAAC' prefix from the "3-AAAC" column sequences if present
guides["3-AAAC"] = guides["3-AAAC"].str.replace('^AAAC', '', regex=True)
print("First few rows of CRISPR guide data after cleaning:")
print(guides.head())

# Load the gRNA frequency data with user input
freq_path = get_input("Enter the path to the gRNA frequency data file: ")
freq_data = pd.read_csv(freq_path, delim_whitespace=True, names=['Frequency', 'Sequence'])
print("First few rows of gRNA frequency data:")
print(freq_data.head())

# Ask user for the low frequency cutoff
low_freq_cutoff = get_input("Enter the low frequency cutoff value: ", int)

# Check sequences from gRNA frequency table against the "3-AAAC" column in guides
matched_data = freq_data[freq_data['Sequence'].isin(guides["3-AAAC"])]

# Merge matched_data with guides to include Gene_name and Guide_Num
merged_data = matched_data.merge(guides[['3-AAAC', 'Gene_name', 'Guide_Num']], left_on='Sequence', right_on='3-AAAC', how='left')

# Calculate the percentage of gRNA sequences that are 100% matched in the guide table
matching_percentage = merged_data['Sequence'].isin(guides["3-AAAC"]).mean() * 100
print(f"Percentage of gRNA sequences that are 100% matched in the CRISPR guide table: {matching_percentage:.2f}%")

# Filter matched sequences to find which have frequencies below the user-defined cutoff
low_freq_matched = merged_data[merged_data['Frequency'] < low_freq_cutoff]

# Print or save the sequences with low frequencies
print(f"Matched guide RNA sequences with frequencies below {low_freq_cutoff}:")
print(low_freq_matched)

# Exporting the low frequency matched sequences to a CSV file
output_csv_path = get_input("Enter the filename to save the low frequency matched sequences: ")
low_freq_matched.to_csv(output_csv_path, index=False)