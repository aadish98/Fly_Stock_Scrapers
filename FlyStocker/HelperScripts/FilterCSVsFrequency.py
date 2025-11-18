import os
import sys
import pandas as pd

def process_csv_files(directory_path):
    # Walk through the directory recursively
    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                try:
                    # Read the CSV into a DataFrame
                    df = pd.read_csv(file_path)

                    # Ensure 'frequency' column exists
                    if 'frequency' in df.columns:
                        # Remove rows where frequency <= 1
                        df = df[df['frequency'] > 1]

                        # Write back to the same file
                        df.to_csv(file_path, index=False)
                        print(f"Processed and updated: {file_path}")
                    else:
                        print(f"'frequency' column not found in: {file_path}")
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Invalid directory path: {directory}")
        sys.exit(1)

    process_csv_files(directory)
