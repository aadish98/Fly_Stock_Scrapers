# Purpose: Filters out FlyBase references with > 50 genes and writes to a blacklist file
# Usage: $python3 LogBlacklistedFlyBaseRefs

import pandas as pd
import glob
import os
import numpy as np


def main():
    # Define the folder path
    folder_path = "../Data/FlyBase/"

    # Get all .tsv files
    tsv_files = glob.glob(os.path.join(folder_path, "*/*.tsv"))

    # Locate the "entity_publication" file
    entity_publication_file = next(
        filter(lambda file: "entity_publication" in file, tsv_files), None
    )
    reference_metadata_file = next(
        filter(lambda file: "fbrf_pmid_pmcid_doi" in file, tsv_files), None
    )
    if not entity_publication_file:
        print("No file containing 'entity_publication' found.")
        return
    if not reference_metadata_file:
        print("No reference_metadata_file found.")
        return

    # Load the entity publication file into a DataFrame
    entity_publication_df = pd.read_csv(
        entity_publication_file, sep="\t", low_memory=False
    )

    reference_metadata_df = pd.read_csv(
        reference_metadata_file, sep="\t", low_memory=False
    )
    print(sorted(np.unique(reference_metadata_df["pub_type"].dropna())))
    # Filter for genes (entity_id starting with FBgn)
    entity_publication_df = entity_publication_df[
        entity_publication_df["entity_id"].str.startswith("FBgn")
    ]
    entity_publication_df = entity_publication_df.drop_duplicates(subset= ["entity_id", "FlyBase_publication_id"])
    # Group by publication ID and count unique entity IDs
    reference_counts = (
        entity_publication_df.groupby("FlyBase_publication_id")["entity_id"]
        .nunique()
        .reset_index(name="unique_gene_count")
    )

    # Split references based on the count threshold
    fewer_than_50 = reference_counts[reference_counts["unique_gene_count"] <= 50]

    more_than_50 = reference_counts[reference_counts["unique_gene_count"] > 50]

    more_than_50 = more_than_50.rename(columns = {"FlyBase_publication_id": "FBrf"})
    more_than_50 = more_than_50.merge(reference_metadata_df, how = "left", on = "FBrf")
    more_than_50 = more_than_50[more_than_50["pub_type"] == "paper"] # Filter on paper type
    more_than_50 = more_than_50.drop_duplicates(subset=["FBrf"])
    more_than_50 = more_than_50.fillna("-")
    more_than_50 = more_than_50.sort_values(by="unique_gene_count", ascending=False)
    # Print results
    print("References with fewer than 50 unique genes:", (fewer_than_50).shape[0], "\n\n")
    print("References with more than 50 unique genes:", (more_than_50).shape[0])
    print("HERERERER")
    print(fewer_than_50.columns)
    print(fewer_than_50[fewer_than_50["FlyBase_publication_id"] == "FBrf0230140"])
    output_path = "../Data/FlyBase/custom/blacklisted_paper_refs.csv"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure the output directory exists
    more_than_50.to_csv(output_path, index=False)
    print("Wrote references with more than 50 genes to ", output_path)
# Run the script
if __name__ == "__main__":
    main()
