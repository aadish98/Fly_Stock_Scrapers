# Example Usage:
# $ python3 GetFBgnIDs.py <InputDirectory> <OutputDirectory> <GeneColumnName (Optional)>
# Handles GO tables (comma separated genes) and simple gene cols
import pandas as pd
import sys
import os
import glob
import numpy as np
import argparse
from itertools import combinations

symbol_to_name = {
    'α': 'alpha',
    'β': 'beta',
    'γ': 'gamma',
    'δ': 'delta',  # Lowercase delta
    'ε': 'epsilon',
    'ζ': 'zeta',
    'η': 'eta',
    'θ': 'theta',
    'ι': 'iota',
    'κ': 'kappa',
    'λ': 'lambda',
    'μ': 'mu',
    'ν': 'nu',
    'ξ': 'xi',
    'ο': 'omicron',
    'π': 'pi',
    'ρ': 'rho',
    'σ': 'sigma',
    'τ': 'tau',
    'υ': 'upsilon',
    'φ': 'phi',
    'χ': 'chi',
    'ψ': 'psi',
    'ω': 'omega',
    'Α': 'Alpha',  # Uppercase Alpha
    'Β': 'Beta',
    'Γ': 'Gamma',
    'Δ': 'Delta',  # Uppercase Delta
    'Ε': 'Epsilon',
    'Ζ': 'Zeta',
    'Η': 'Eta',
    'Θ': 'Theta',
    'Ι': 'Iota',
    'Κ': 'Kappa',
    'Λ': 'Lambda',
    'Μ': 'Mu',
    'Ν': 'Nu',
    'Ξ': 'Xi',
    'Ο': 'Omicron',
    'Π': 'Pi',
    'Ρ': 'Rho',
    'Σ': 'Sigma',
    'Τ': 'Tau',
    'Υ': 'Upsilon',
    'Φ': 'Phi',
    'Χ': 'Chi',
    'Ψ': 'Psi',
    'Ω': 'Omega'
}

def replace_symbol(gene_series):
    for symbol, name in symbol_to_name.items():
        gene_series = gene_series.str.replace(symbol, name, regex=False)
    return gene_series

def create_expanded_mappings(mappings_df):
    """
    Create expanded mappings for each relevant column to ensure all entries are processed
    (e.g., pipe-separated values in 'fullname_synonym(s)' and 'symbol_synonym(s)').
    """
    def expand_synonyms(column_name):
        """
        Expand a column containing pipe-separated synonyms into individual mappings.
        """
        if column_name not in mappings_df.columns:
            raise ValueError(f"Column '{column_name}' not found in the DataFrame.")

        # Select relevant rows where the column is not null
        relevant_rows = mappings_df.dropna(subset=[column_name]).copy()

        # Split pipe-separated values into lists
        relevant_rows[column_name] = relevant_rows[column_name].str.split('|')

        # Explode to create individual rows for each synonym
        expanded = relevant_rows.explode(column_name).rename(columns={column_name: 'synonym'})
        expanded['synonym'] = expanded['synonym'].str.strip()  # Remove whitespace

        return expanded[['synonym', 'primary_FBid']]

    # Create mappings for each relevant column
    columns_to_expand = ['current_symbol', 'current_fullname', 'fullname_synonym(s)', 'symbol_synonym(s)']
    all_mappings = pd.DataFrame()

    for column in columns_to_expand:
        print(f"Processing column: {column}")
        expanded_mapping = expand_synonyms(column)
        all_mappings = pd.concat([all_mappings, expanded_mapping], ignore_index=True)

    # Drop duplicates to ensure unique mappings
    all_mappings = all_mappings.drop_duplicates()

    # Convert to a dictionary for efficient lookups
    synonym_to_fbgnid_map = dict(zip(all_mappings['synonym'], all_mappings['primary_FBid']))

    return synonym_to_fbgnid_map


def map_gene_ids(df, gene_to_fbgnid_main, gene_to_fbgnid_synonym, gene_col):
    """
    Efficiently maps ext_gene in df to flybase_gene_id using main and synonym mappings,
    while applying prioritized step combinations in a vectorized manner.
    """
    import re

    def clean_gene_vectorized(series, steps):
        """
        Apply a sequence of cleaning steps to a pandas Series in a vectorized manner.
        """
        for step in steps:
            if step == 'lowercase':
                series = series.str.lower()
            elif step == 'remove_zeros':
                series = series.str.rstrip('0').str.replace('-0', '-', regex=False).str.replace('-00', '-', regex=False)
            elif step == 'capitalize':
                series = series.str.capitalize()
            elif step == 'remove_hyphen':
                series = series.str.replace('-', '', regex=False)
            elif step == 'number_to_end':
                series = series.str.replace(r"^(\d+)-(\D+)$", r"\2\1", regex=True)  # Handles cases like "1-Sep" -> "Sep1"
            elif step == 'number_to_start':
                series = series.str.replace(r"^(\D+)-(\d+)$", r"\2\1", regex=True)  # Handles cases like "Sep-1" -> "1Sep"
            elif step == 'add_cr':
                series = "CR" + series
            elif step == 'uppercase':
                series = series.str.upper()
        return series

    # Step 1: Initial mapping
    df[gene_col] = df[gene_col].str.strip()
    df['flybase_gene_id'] = df[gene_col].map(gene_to_fbgnid_main)

    # Step 2: Synonym mapping for initially unmapped genes
    unmapped_genes_mask = df['flybase_gene_id'].isna()
    print(f'TOTAL GENES: {df.shape[0]}')
    print(f'Initially missed genes: {unmapped_genes_mask.sum()}')
    df.loc[unmapped_genes_mask, 'flybase_gene_id'] = df.loc[unmapped_genes_mask, gene_col].map(gene_to_fbgnid_synonym)
    print(f'Unmatched synonyms: {df["flybase_gene_id"].isna().sum()}')
    print(f'{df[df["flybase_gene_id"].isna()]}')
    # Step 3: Prioritized step combinations
    steps = ['lowercase', 'remove_zeros', 'remove_hyphen', 'number_to_end', 'uppercase'] # 'add_cr'
    all_step_combinations = [combo for r in range(1, len(steps) + 1) for combo in combinations(steps, r)]
    # for r in range(1, len(steps) + 1):
    #     all_step_combinations.extend(permutations(steps, r))

    # Cache intermediate results
    cache = {}
    prev_unmapped_count = df['flybase_gene_id'].isna().sum()

    # Step 4: Apply combinations iteratively
    for step_combo in all_step_combinations:
        unmapped_genes_mask = df['flybase_gene_id'].isna()
        if not unmapped_genes_mask.any():
            break  # Exit if all genes are mapped

        # Process only unmapped genes
        unmapped_genes = df.loc[unmapped_genes_mask, gene_col]
        if step_combo in cache:
            cleaned_genes = cache[step_combo]
        else:
            cleaned_genes = clean_gene_vectorized(unmapped_genes, step_combo)
            cache[step_combo] = cleaned_genes

        # Try mapping with main and synonym mappings
        df.loc[unmapped_genes_mask, 'flybase_gene_id'] = cleaned_genes.map(gene_to_fbgnid_main)
        unmapped_genes_mask = df['flybase_gene_id'].isna()
        df.loc[unmapped_genes_mask, 'flybase_gene_id'] = cleaned_genes.map(gene_to_fbgnid_synonym)

        current_unmapped_count = df['flybase_gene_id'].isna().sum()
        if current_unmapped_count == prev_unmapped_count:
            continue  # Skip redundant combinations
        else:
            print(f'{step_combo} was succesful in resolving {prev_unmapped_count - current_unmapped_count} genes')
        prev_unmapped_count = current_unmapped_count

    # Step 5: Log remaining unmapped genes
    remaining_unmapped_mask = df['flybase_gene_id'].isna()
    print(f"Remaining unmapped genes after all steps: {remaining_unmapped_mask.sum()}")
    if remaining_unmapped_mask.sum() > 0:
        print("Remaining unmapped genes:")
        print(df.loc[remaining_unmapped_mask, gene_col])

    return df



# Main script
if __name__ == "__main__":
    # Get directory paths and optional gene column name from command-line arguments
    print('\n')
    print("HERE NEW RUN NEW RUN NEW RUN ===============================================================")
    directory_path = sys.argv[1]  # Directory containing the input CSV files
    default_gene_col = "ext_gene"
    gene_column = sys.argv[2] if len(sys.argv) > 2 else default_gene_col  # Column with gene names, default is 'ext_gene'
    print("HERE")
    print(gene_column)
    csv_files = glob.glob(os.path.join(directory_path, '**.csv'), recursive=True)
    excel_files = glob.glob(os.path.join(directory_path, '**.xlsx'), recursive=True)
    # Load the FlyBase synonyms mapping file
    columns_as_strings = {
    'current_symbol': str,
    'current_fullname': str,
    'fullname_synonym(s)': str,
    'symbol_synonym(s)': str,
    'primary_FBid': str
    }


    mappings_df = pd.read_csv('~/Desktop/Projects.nosync/Allada-Lab/Fly_Stock_Scrapers/Data/FlyBase/Genes/fb_synonym_fb_2025_03.tsv', sep='\t', header=0, low_memory=False, dtype=columns_as_strings)
    mappings_df = mappings_df[mappings_df.organism_abbreviation == "Dmel"]
    # Keep only entries where the primary gene ID starts with 'FBgn'
    mappings_df = mappings_df[mappings_df['primary_FBid'].str.startswith('FBgn', na=False)]
    synonym_to_fbgnid_map = create_expanded_mappings(mappings_df)

    gene_to_fbgnid_main = dict(zip(mappings_df['current_symbol'].str.strip(), mappings_df['primary_FBid']))

    # Process each CSV file
    for csv_file in csv_files:
        if len(csv_files) == 0:
            break
        print(f"Processing file: {csv_file}")

        # Load the CSV file into a DataFrame
        df = pd.read_csv(csv_file, encoding = 'utf-8-sig', na_values=[""], dtype={gene_column: str}, low_memory=False)

        if gene_column not in df.columns:
            print(f"'{gene_column}' column not found in {csv_file}. Skipping file.")
            continue  # Skip to the next file if the gene column is not found

        # Remove 'flybase_gene_id' column if it exists (to avoid conflicts)
        if 'flybase_gene_id' in df.columns:
            df.drop(columns=['flybase_gene_id'], inplace=True)
        if ("GO" in csv_file and "Genes" in df.columns):
            df['Genes'] = df['Genes'].str.split(', ')
            # Explode the DataFrame so each gene gets its own row
            df = df.explode('Genes').reset_index(drop=True)
        df[gene_column + "_new"] = replace_symbol(df[gene_column])
        # Drop empty columns
        df = df.dropna(subset=[gene_column + "_new"])
        # Map the genes in the specified column to FlyBase gene IDs
        df = map_gene_ids(df, gene_to_fbgnid_main, synonym_to_fbgnid_map, gene_column+ "_new")
        # Save the resulting DataFrame to the output directory
        df.fillna("-", inplace = True)
        df.drop(columns=[gene_column + "_new"], inplace= True)
        df.to_csv(csv_file, index=False, encoding = 'utf-8-sig')

        print(f"File saved to {csv_file} with {df.shape[0]} rows.")

    print('All files processed.')
