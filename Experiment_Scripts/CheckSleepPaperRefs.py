import pandas as pd
import sys
import os
import glob


# Function to check references and map to the new column
def get_sleep_references(row, mappings_set):
    # Split the references in the PAPER_REFERENCES column
    references = row.split("; ")

    # Filter only the references that are found in the mappings set
    matched_refs = [ref for ref in references if ref in mappings_set]

    # Return the matched references, or 'NONE' if no matches found
    return "; ".join(matched_refs) if matched_refs else "NONE"


# Main script
if __name__ == "__main__":
    # Get directory paths from command-line arguments
    directory_path = sys.argv[1]  # Directory containing the input CSV files

    # Load the mapping DataFrame and convert 'reference_id' to a set for fast lookup
    mappings_df = pd.read_csv(
        "Data/FlyBase/QueryResults/Sleep/TitleAbstract_Sleep_PubType_Paper_RefIDs.csv"
    )
    mappings_set = set(
        mappings_df["reference_id"].tolist()
    )  # Optimized for fast lookups

    # Load all CSV files from the directory
    csv_files = [
        f
        for f in glob.glob(os.path.join(directory_path, "**.csv"), recursive=True)
        if not os.path.basename(f).startswith("~$")
    ]

    # Process each CSV file
    for csv_file in csv_files:
        print(f"Processing file: {csv_file}")

        # Open the file in chunks
        chunk_size = 1000  # Adjust chunk size based on your memory capacity
        chunks = []

        # Process the file chunk by chunk
        for chunk in pd.read_csv(csv_file, chunksize=chunk_size, low_memory=False):
            # Apply the get_sleep_references function (optimized using the set)
            chunk["SLEEP_PAPER_REFS"] = chunk["PAPER_REFERENCES"].apply(
                lambda x: get_sleep_references(x, mappings_set)
            )
            chunks.append(chunk)

        # Concatenate all processed chunks and write back to CSV
        df_concat = pd.concat(chunks)
        df_concat.to_csv(csv_file, index=False)

        # Output the number of matching references
        num_matching_stocks = df_concat[df_concat["SLEEP_PAPER_REFS"] != "NONE"].shape[
            0
        ]
        print(
            f"File saved to {csv_file} with {df_concat.shape[0]} stocks, of which {num_matching_stocks} stocks are in sleep papers\n\n"
        )

    print("All files processed.")
