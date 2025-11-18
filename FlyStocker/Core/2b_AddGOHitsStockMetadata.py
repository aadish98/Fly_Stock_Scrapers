import pandas as pd
import os
import glob
import sys
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import argparse

from AddStockMetadata import get_sleep_genes, get_neurodegen_genes, get_sleep_paper_refs, get_nd_paper_refs
# Ge
# Function to remove characters up to and including '~', then sort terms in each row
def clean_and_sort_term(term):
    # Split by commas, remove up to '~', sort alphabetically, and join back
    sorted_terms = sorted(set([t.split('~')[-1].strip() for t in term.split(',')]))
    return ','.join(sorted_terms)

# Function to extract all unique GO terms across the dataset
def get_unique_go_terms(df, term_col='GO_Terms'):
    all_terms = set()
    for terms in df[term_col]:
        term_list = terms.replace('\n', ',').split(',')
        all_terms.update([t.strip() for t in term_list])
    return sorted(all_terms)

# Function to create one-hot encoding for GO terms
def create_one_hot_encoding(df, unique_terms, term_col='GO_Terms'):
    for term in unique_terms:
        df[term] = df[term_col].apply(lambda x: 1 if term in x else 0)
    return df

if __name__ == "__main__":
    # Get directory paths and optional gene column name from command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("directory_path")
    parser.add_argument("target_dir")
    parser.add_argument("-neurodegen", action="store_true", help="Set to true if neurodegeneration is required")
    parser.add_argument("-sleep", action="store_true", help="Set to true if neurodegeneration is required")
    parser.add_argument("-sortby", choices=["frequency"], help="Specify sorting criteria")
    args = parser.parse_args()

    directory_path = args.directory_path
    target_dir = args.target_dir
    neurodegen_flag = args.neurodegen
    sleep_flag = args.sleep
    sortby = args.sortby
    csv_files = glob.glob(os.path.join(directory_path, '**.csv'), recursive=True)
    df_list = []

    # Loop through the list of CSV files and read each into a DataFrame
    for file in csv_files:
        df = pd.read_csv(file, keep_default_na=False)  # Read each CSV file
        df_list.append(df)  # Append the DataFrame to the list

    # Concatenate all the DataFrames in the list into a single DataFrame
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df['Genes'] = combined_df['Genes'].str.split(',').apply(lambda x: [gene.strip() for gene in x])
    combined_df = combined_df.explode('Genes', ignore_index=True)
    combined_df = combined_df.groupby('Genes', as_index=False).agg({
        'Term': lambda x: ', '.join(sorted(set(x)))  # Combine unique terms, sorted, into a single string
    })
    # Rename the 'Term' column to 'GO_Term'
    combined_df.rename(columns={'Term': 'GO_Terms'}, inplace=True)
    if ("Genes" not in combined_df.columns):
        raise RuntimeError("METADATA DF DOESNT CONTAIN FBgnID!!!")
    target_files = glob.glob(os.path.join(target_dir, "**.csv"), recursive=True)

    for target_file in target_files:
        print(target_file)
        target_df = pd.read_csv(target_file, keep_default_na=False)
        print(f'Original file length: {target_df.shape[0]}')
        # Add flybase gene id for stock table

        if "GeneID_x" in target_df.columns:
            target_df["GeneID"] = target_df["GeneID_x"]
            target_df.drop(columns=["GeneID_x", "GeneID_y"])
        if "flybase_gene_id_x" in target_df.columns:
            target_df["flybase_gene_id"] = target_df["flybase_gene_id_x"]
            target_df.drop(columns=["flybase_gene_id_x", "flybase_gene_id_y"])
        if "GO_Terms_x" in target_df.columns:
            target_df.drop(columns=["GO_Terms_x"], inplace = True)
        if "GO_Terms_y" in target_df.columns:
            target_df.drop(columns=["GO_Terms_y"], inplace = True)
        if "GO_Terms" in target_df.columns:
            target_df.drop(columns=["GO_Terms"], inplace = True)

        original_columns = target_df.columns.to_list()
        meta_df = pd.merge(target_df, combined_df[["Genes", "GO_Terms"]], how="left", left_on="flybase_gene_id", right_on="Genes")
        if "flybase_gene_id_x" in meta_df.columns:
            meta_df["flybase_gene_id"] = meta_df["flybase_gene_id_x"]
            meta_df.drop(columns=["flybase_gene_id_x", "flybase_gene_id_y"])
        if "FDR_x" in meta_df.columns:
            meta_df["FDR"] = meta_df["FDR_x"]
            meta_df.drop(columns=["FDR_x", "FDR_y"])

        # Add genes relevant publications
        if sleep_flag:
            meta_df["GENE_HAS_SLEEP_PAPER_REF"] = meta_df["flybase_gene_id"].apply(
                get_sleep_genes
            )
            meta_df["GENE_SLEEP_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                get_sleep_paper_refs
            )
        if neurodegen_flag:
            meta_df["GENE_HAS_NEURODEGEN_PAPER_REF"] = meta_df["flybase_gene_id"].apply(
                get_neurodegen_genes
            )
            meta_df["GENE_NEURODEGEN_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                get_nd_paper_refs
            )

        # Reorder the columns
        cols_to_order = [col for col in ["GENE_HAS_PAPER_REFS", "GENE_HAS_SLEEP_PAPER_REF", "GENE_HAS_NEURODEGEN_PAPER_REF", "GO_Terms"] if col in meta_df.columns]
        other_cols = [col for col in meta_df.columns if col not in cols_to_order]
        meta_df = meta_df[other_cols + cols_to_order]
        meta_df.drop(columns=['Genes'], inplace=True)

        meta_df['GO_Terms'] = meta_df['GO_Terms'].fillna('-')
        meta_df.to_csv(target_file, index=False)
        print(f'Final length: {meta_df.shape[0]}')
