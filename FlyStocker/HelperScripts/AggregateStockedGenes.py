import os
import glob
import pandas as pd


def aggregate_csv_files(input_dir, output_filename):
    """
    Recursively finds all CSV files in the directory (excluding those with 'summary' in their name),
    extracts columns ASSOCIATED_GENE, flybase_gene_id, and StockID, aggregates them,
    and saves the result in the same directory.
    """
    csv_files = glob.glob(os.path.join(input_dir, '**/*.csv'), recursive=True)

    # Filter out files containing 'summary' in the filename
    csv_files = [f for f in csv_files if 'summary' not in os.path.basename(f).lower()]

    if not csv_files:
        print("No matching CSV files found.")
        return

    data_frames = []

    for file in csv_files:
        try:
            df = pd.read_csv(file, usecols=['ASSOCIATED_GENE', 'flybase_gene_id', 'StockID'], dtype=str,
                             low_memory=False)
            data_frames.append(df)
            print(f"Processed: {file}")
        except Exception as e:
            print(f"Skipping {file} due to error: {e}")

    if not data_frames:
        print("No valid data to aggregate.")
        return

    # Concatenate all DataFrames
    aggregated_df = pd.concat(data_frames, ignore_index=True).drop_duplicates()

    # Save the aggregated DataFrame
    output_path = os.path.join(input_dir, output_filename)
    aggregated_df.to_csv(output_path, index=False, encoding='utf-8-sig')
    print(f"Aggregated data saved to: {output_path}")


if __name__ == "__main__":
    input_directory = "/Users/aadishms/University of Michigan Dropbox/Aadish Shah/Computational_group/Aadish's Brain Seq Analysis/FlyBase Stocks/TDP43"
    output_file = "TDP43_Aggregated.csv"
    aggregate_csv_files(input_directory, output_file)