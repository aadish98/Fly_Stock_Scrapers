import pandas as pd
import os
import glob
import argparse


def create_fbgn_to_extgene_map(mappings_df):
    """
    Create a dictionary mapping FlyBase gene IDs (FBgn) to primary gene symbols.
    """
    if 'primary_FBid' not in mappings_df.columns or 'current_symbol' not in mappings_df.columns:
        raise ValueError("Mapping file must contain 'primary_FBid' and 'current_symbol' columns.")
    return dict(zip(mappings_df['primary_FBid'], mappings_df['current_symbol'].str.strip()))


def update_missing_ext_gene(df, fbgn_to_extgene_map, fbgn_col='flybase_gene_id', ext_gene_col='ext_gene'):
    """
    Update the ext_gene column only for rows where ext_gene is '-' using the FlyBase mapping.
    """
    # Identify rows where ext_gene equals '-'
    missing_mask = df[ext_gene_col] == '-'

    if missing_mask.any():
        # For rows where ext_gene is '-', map the FBgnID to its corresponding gene symbol.
        df.loc[missing_mask, ext_gene_col] = df.loc[missing_mask, fbgn_col].map(
            lambda fbgn: fbgn_to_extgene_map.get(fbgn, '-'))
        print(f"Updated {missing_mask.sum()} rows in '{ext_gene_col}' column where value was '-'.")
    else:
        print(f"No rows with '-' in '{ext_gene_col}' found. No update needed.")

    return df


def process_excel_files(input_dir, fbgn_col = 'flybase_gene_id'):
    """
    Process Excel files in the given directory and update the 'ext_gene' column
    in the 'aggregated_results' sheet only for rows where the value is '-'.
    """
    # Find Excel files (.xlsx and .xls) recursively
    excel_files = glob.glob(os.path.join(input_dir, '**/*.xlsx'), recursive=True) + \
                  glob.glob(os.path.join(input_dir, '**/*.xls'), recursive=True)

    if not excel_files:
        print("No Excel files found in the specified directory.")
        return

    # Load the FlyBase mappings file
    mapping_file = os.path.expanduser(
        '~/Desktop/Projects.nosync/Allada-Lab/Fly_Stock_Scrapers/Data/FlyBase/Genes/fb_synonym_fb_2024_06.tsv')
    mappings_df = pd.read_csv(mapping_file, sep='\t', header=0, dtype={'primary_FBid': str, 'current_symbol': str})
    fbgn_to_extgene_map = create_fbgn_to_extgene_map(mappings_df)

    # Process each Excel file
    for excel_file in excel_files:
        print(f"\nProcessing file: {excel_file}")
        try:
            sheets = pd.read_excel(excel_file, sheet_name=None)
        except Exception as e:
            print(f"Error reading file {excel_file}: {e}")
            continue

        if 'aggregated_results' not in sheets:
            print(f"Sheet 'aggregated_results' not found in {excel_file}. Skipping this file.")
            continue

        df = sheets['aggregated_results']

        # Check for the required columns
        if fbgn_col not in df.columns:
            print(
                f"Column '{fbgn_col}' not found in sheet 'aggregated_results' in {excel_file}. Skipping update for this sheet.")
            continue
        if 'ext_gene' not in df.columns:
            print(
                f"Column 'ext_gene' not found in sheet 'aggregated_results' in {excel_file}. Skipping update for this sheet.")
            continue

        # Update only the rows where ext_gene is '-'
        df = update_missing_ext_gene(df, fbgn_to_extgene_map, fbgn_col=fbgn_col, ext_gene_col='ext_gene')
        sheets['aggregated_results'] = df

        # Write all sheets back to the Excel file
        try:
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                for sheet_name, data in sheets.items():
                    data.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"Updated file saved: {excel_file}")
        except Exception as e:
            print(f"Error saving file {excel_file}: {e}")

    print("\nProcessing complete.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Update the ext_gene column in Excel files (aggregated_results sheet) only for rows with '-' using FlyBase mapping."
    )
    parser.add_argument("input_dir", help="Directory containing Excel files")
    parser.add_argument("--gene_col", default="flybase_gene_id",
                        help="Name of the column with FlyBase gene IDs (default: flybase_gene_id)")
    args = parser.parse_args()
    process_excel_files(args.input_dir, args.gene_col)
