import os
import sys
import pandas as pd
import re

# Function to extract the gene name after "RNAi of ..."
def extract_predicted_gene(comment):
    if pd.isna(comment):
        return None  # Return None if the comment is missing
    match = re.search(r"RNAi of ([^\s]+)", comment)  # Capture everything after "RNAi of " until the next space
    return match.group(1) if match else None

# Ensure command-line argument is provided
if len(sys.argv) < 2:
    print("Usage: python map_stocks.py <directory_path>")
    sys.exit(1)

# Get the directory from command-line argument
input_dir = sys.argv[1]

# Verify directory exists
if not os.path.isdir(input_dir):
    print(f"Error: Directory '{input_dir}' does not exist.")
    sys.exit(1)

# Path to the BDSC data file
BDSC_FILE = "/Users/aadishms/Desktop/Projects.nosync/Allada-Lab/Fly_Stock_Scrapers/Data/BDSC/stockcomps_map_comments.csv"

# Read the BDSC data file
bdsc_df = pd.read_csv(
    BDSC_FILE, dtype={"Stk #": str}, usecols=["Stk #", "comment1", "comment2", "comment3"], low_memory=False
)

# Extract predicted gene
bdsc_df["ext_gene"] = bdsc_df["comment1"].apply(extract_predicted_gene)

# Create the output directory
output_dir = os.path.join(input_dir, "Mapped_Stocks")
os.makedirs(output_dir, exist_ok=True)

# Find all `.xlsx` and `.csv` files in the input directory (excluding temp Excel files)
files_to_process = [f for f in os.listdir(input_dir) if f.endswith(('.xlsx', '.csv')) and not f.startswith('~$')]

for file in files_to_process:
    input_file_path = os.path.join(input_dir, file)
    print(f"Processing {file}...")

    try:
        # Determine file type and read accordingly
        if file.endswith(".xlsx"):
            df = pd.read_excel(input_file_path, dtype={"StockID": str}, engine='openpyxl')
        else:  # file is .csv
            df = pd.read_csv(input_file_path, dtype={"StockID": str})
        df = df.dropna(subset=['StockID'])
        # DONT MAP CONTROL LINES
        df = df[~df['StockID'].isin(['36303', '36304'])]
        df = df.dropna(subset=['StockID'])
        #df['StockID'] = pd.to_numeric(df['StockID'], errors='coerce')

        df["StockID"] = df.StockID.astype(str)
        # Merge with the BDSC data
        merged_df = df.merge(bdsc_df, left_on="StockID", right_on="Stk #", how="left")

        # Drop the 'Stk #' column after merging
        columns_to_drop = ['Stk #', 'comment2', 'comment3'] + [col for col in merged_df.columns if 'Unnamed' in col]
        merged_df.drop(columns=columns_to_drop, inplace=True, errors='ignore')

        # Construct output file path
        output_file_path = os.path.join(output_dir, file)

        # Save in the same format as the input
        if file.endswith(".xlsx"):
            merged_df.to_excel(output_file_path, index=False, engine='openpyxl')
        else:
            merged_df.to_csv(output_file_path, index=False)

        print(f"Saved mapped file as {output_file_path}")

    except Exception as e:
        print(f"Error processing {file}: {e}")

print("Processing complete!")
