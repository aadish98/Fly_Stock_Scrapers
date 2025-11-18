import pandas as pd
import os
import glob
import sys
import argparse
# 1 cmdline arg: directory_path

def createAlleleFolderPath(directory_path):
    parent_directory = os.path.dirname(os.path.join(directory_path))
    print(parent_directory)
    # Define the name of the new subfolder to create
    new_folder_name = 'AlleleIDs'

    # Create the full path for the new folder
    new_folder_path = os.path.join(parent_directory, new_folder_name)

    # Create the new folder if it doesn't already exist
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        print(f"Created new folder: {new_folder_path}")

    return new_folder_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory_path")
    args = parser.parse_args()
    directory_path = args.directory_path
    allele_map = "~/Desktop/Projects.nosync/Allada-Lab/Fly_Stock_Scrapers/Data/FlyBase/Alleles_And_Stocks/fbal_to_fbgn_fb_2024_05.tsv"
    allele_map_df = pd.read_csv(allele_map, delimiter='\t', low_memory=False)
    csv_files = glob.glob(os.path.join(directory_path, '**.csv'), recursive=True)
    alleleID_folder_path = createAlleleFolderPath(directory_path)
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        merged_df = df.merge(allele_map_df, how = 'left', left_on= 'flybase_gene_id', right_on = 'GeneID')
        alleleID_file_name = os.path.basename(csv_file)
        merged_df.drop(columns = ["GeneID",	"GeneSymbol"], inplace= True)
        merged_df.drop_duplicates(subset = ["AlleleID"], inplace= True)
        # Split into chunks of 1000 rows
        num_chunks = (len(merged_df) - 1) // 1000 + 1
        for i in range(num_chunks):
            chunk = merged_df.iloc[i * 1000:(i + 1) * 1000]
            chunk_file_name = alleleID_file_name.replace(".csv", f"_part{i + 1}.csv")
            chunk.to_csv(os.path.join(alleleID_folder_path, chunk_file_name), index=False)
        print(f'Saved {num_chunks} chunk(s) for file {alleleID_file_name} ({merged_df.shape[0]} alleles)')

