import pandas as pd
import sys
import os
import glob
import argparse


def create_fbgn_to_extgene_map(mappings_df):
    """
    Create a dictionary mapping FlyBase gene IDs (FBgn) to external gene symbols (ext_gene).
    """
    if 'primary_FBid' not in mappings_df.columns or 'current_symbol' not in mappings_df.columns:
        raise ValueError("Missing required columns in mappings file.")

    fbgn_to_extgene_map = dict(zip(mappings_df['primary_FBid'], mappings_df['current_symbol'].str.strip()))
    return fbgn_to_extgene_map


def map_fbgn_to_extgene(df, fbgn_to_extgene_map, gene_col):
    """
    Map flybase_gene_id to ext_gene using the provided dictionary.
    If multiple FBgn IDs exist in a single row, map each separately and join the results.
    """
    df[gene_col] = df[gene_col].astype(str).str.strip()

    def map_gene_list(fbgn_list):
        """Convert FBgnID list to gene symbols while handling missing mappings."""
        genes = [fbgn_to_extgene_map.get(fbgn.strip(), fbgn.strip()) for fbgn in fbgn_list.split(", ")]
        return ", ".join(genes)  # Return as a comma-separated string

    df['Genes'] = df[gene_col].apply(map_gene_list)

    # Log unmatched genes
    unmatched_genes = df[df['Genes'].isna()][gene_col]
    print(f"Unmatched genes count: {len(unmatched_genes)}")
    if len(unmatched_genes) > 0:
        print("Unmatched genes:")
        print(unmatched_genes.unique())

    return df


def process_files(input_dir, gene_column):
    """
    Process CSV files in the given directory and convert FlyBase gene IDs to external gene symbols.
    """
    csv_files = glob.glob(os.path.join(input_dir, '**/*.csv'), recursive=True)

    if not csv_files:
        print("No CSV files found in the directory.")
        return

    # Load FlyBase mappings file
    mappings_df = pd.read_csv(
        '~/Desktop/Projects.nosync/Allada-Lab/Fly_Stock_Scrapers/Data/FlyBase/Genes/fb_synonym_fb_2024_06.tsv',
        sep='\t', header=0, low_memory=False, dtype={'primary_FBid': str, 'current_symbol': str})
    fbgn_to_extgene_map = create_fbgn_to_extgene_map(mappings_df)

    for csv_file in csv_files:
        print(f"Processing file: {csv_file}")
        df = pd.read_csv(csv_file, encoding='utf-8-sig', dtype={gene_column: str})

        if gene_column not in df.columns:
            print(f"'{gene_column}' column not found in {csv_file}. Skipping.")
            continue

        # Remove 'Genes' column if it exists (avoid conflicts)
        if 'Genes' in df.columns:
            df.drop(columns=['Genes'], inplace=True)

        df = map_fbgn_to_extgene(df, fbgn_to_extgene_map, gene_column)

        # Save updated DataFrame
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"Updated file saved: {csv_file} with {df.shape[0]} rows.")

    print("All files processed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FlyBase gene IDs to external gene symbols.")
    parser.add_argument("input_dir", help="Directory containing input CSV files.")
    parser.add_argument("--gene_col", default="flybase_gene_id",
                        help="Column containing FlyBase gene IDs (default: flybase_gene_id)")

    args = parser.parse_args()
    process_files(args.input_dir, args.gene_col)
