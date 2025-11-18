# Usage: python3 3_SummarizeProcessedFlybaseTsv.py /path/to/process_csv_files /path/to/config.json
import pandas as pd
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import warnings
import statsmodels.api as sm
from statsmodels.graphics.mosaicplot import mosaic
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from adjustText import adjust_text
import yaml
from itertools import product
import json
import argparse
import re
import csv
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from HelperScripts.plotting_functions import *
from HelperScripts.processing_helper_functions import *
import time

def write_to_excel_autosized(file_path, dataframes, sheet_names):
    """
    Write multiple DataFrames to an Excel file with autosized column widths.

    Parameters:
    - file_path (str): Path to save the Excel file.
    - dataframes (list of pd.DataFrame): List of DataFrames to write to Excel.
    - sheet_names (list of str): List of sheet names corresponding to each DataFrame.
    """
    if len(dataframes) != len(sheet_names):
        raise ValueError("Number of DataFrames and sheet names must match.")

    with pd.ExcelWriter(file_path, engine="xlsxwriter") as writer:
        for df, sheet_name in zip(dataframes, sheet_names):
            # Write the DataFrame to the Excel sheet
            df.to_excel(writer, sheet_name=sheet_name, index=False)

            # Access the worksheet and workbook
            workbook = writer.book
            worksheet = writer.sheets[sheet_name]

            # Autosize columns
            for i, column in enumerate(df.columns):
                max_length = max(
                    df[column].astype(str).map(len).max(),  # Length of the data
                    len(column),  # Length of the column header
                )
                worksheet.set_column(
                    i, i, max_length + 2
                )  # Add padding for readability

    print(f"Excel file saved with autosized columns at: {file_path}")


def read_insertion_mapping_df():
    insertion_mapping_df = pd.read_csv(
        "../../Data/FlyBase/Transgenic_Contructs_And_Insertions/insertion_mapping_fb_2024_05.tsv",
        sep="\t",
    )
    return insertion_mapping_df


def read_gene_mapping_df():
    gene_mapping_df = pd.read_csv(
        "../../Data/FlyBase/Genes/gene_map_table_fb_2024_05.tsv", sep="\t"
    )
    gene_mapping_df = gene_mapping_df[
        gene_mapping_df["organism_abbreviation"] == "Dmel"
    ]
    gene_mapping_df = gene_mapping_df.astype(str)
    return gene_mapping_df


def limit_rows_per_gene(df, stocksPerGene):
    return (
        df.groupby("ASSOCIATED_GENE")
        .apply(lambda x: x.sample(n=min(len(x), stocksPerGene)))
        .reset_index(drop=True)
    )


def affected_chromosomes(bloom_stock_no, bdsc_stocks):
    result = bdsc_stocks[bdsc_stocks["Stk #"] == bloom_stock_no][
        "All Affected Chromosomes"
    ]
    return result.iloc[0] if not result.empty else ""


def extract_affected_genes(genotype):
    if not isinstance(genotype, str):
        return "-"
    # Split genotype by delimiters and check each token against the set of gene symbols
    affected_genes = [gene for gene in all_gene_symbols if str(gene) in genotype]
    # Join the matches into a comma-separated string
    return ", ".join(map(str, affected_genes)) if affected_genes else "-"


def load_split_config(config_path):
    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def apply_filters(df, filters_config):
    filters = {}
    for filter_name, condition in filters_config.items():
        if condition["type"] == "contains":
            filters[filter_name] = df[condition["column"]].str.contains(
                condition["value"], case=True, na=False
            )
        elif condition["type"] == "doesnt_contain":
            filters[filter_name] = ~df[condition["column"]].str.contains(
                condition["value"], case=True, na=False
            )
        elif condition["type"] == "equals":
            filters[filter_name] = df[condition["column"]] == condition["value"]
        elif condition["type"] == "doesnt_equal":
            filters[filter_name] = df[condition["column"]] != condition["value"]
        elif condition["type"] == "insertion_site":
             filters[filter_name] = (
                df[condition["column"]].str.contains(condition["value"], case=True, na=False)
                & (df["multiple_insertions"] == False)
            )
        # Add more conditions as needed
    return filters


def find_overlaps(stock_refs, gene_refs):
    if stock_refs == "-" or gene_refs == "-":
        return "-"
    
    # Split and strip references
    stock_set = set(ref.strip() for ref in stock_refs.split(", "))
    gene_set = set(ref.strip() for ref in gene_refs.split(", "))
    
    # Find overlap
    overlap = stock_set.intersection(gene_set)
    return ", ".join(sorted(overlap)) if overlap else "-"


def calculate_relevance_score(stock_sleep_papers, stock_paper_refs):
    if (
        stock_sleep_papers != "-"
    ):  # Allele sleep papers exist
        return 2
    elif (
        stock_paper_refs != "-"
    ):  # No overlaps, but STOCK_PAPER_REFS is not empty
        return 1
    else:  # STOCK_PAPER_REFS is empty
        return 0

def write_to_excel_autosized(file_path, dataframes, sheet_names):
    """
    Write multiple DataFrames to an Excel file with autosized column widths.

    Parameters:
    - file_path (str): Path to save the Excel file.
    - dataframes (list of pd.DataFrame): List of DataFrames to write to Excel.
    - sheet_names (list of str): List of sheet names corresponding to each DataFrame.
    """
    if len(dataframes) != len(sheet_names):
        raise ValueError("Number of DataFrames and sheet names must match.")

    with pd.ExcelWriter(file_path, engine="xlsxwriter") as writer:
        for df, sheet_name in zip(dataframes, sheet_names):
            # Write the DataFrame to the Excel sheet
            df.to_excel(writer, sheet_name=sheet_name, index=False)

            # Access the worksheet and workbook
            workbook = writer.book
            worksheet = writer.sheets[sheet_name]

            # Autosize columns
            for i, column in enumerate(df.columns):
                max_length = max(
                    df[column].astype(str).map(len).max(),  # Length of the data
                    len(column),  # Length of the column header
                )
                worksheet.set_column(
                    i, i, max_length + 2
                )  # Add padding for readability

    print(f"Excel file saved with autosized columns at: {file_path}")


def read_insertion_mapping_df():
    insertion_mapping_df = pd.read_csv(
        "../../Data/FlyBase/Transgenic_Contructs_And_Insertions/insertion_mapping_fb_2024_05.tsv",
        sep="\t",
    )
    return insertion_mapping_df


def read_gene_mapping_df():
    gene_mapping_df = pd.read_csv(
        "../../Data/FlyBase/Genes/gene_map_table_fb_2024_05.tsv", sep="\t"
    )
    gene_mapping_df = gene_mapping_df[
        gene_mapping_df["organism_abbreviation"] == "Dmel"
    ]
    gene_mapping_df = gene_mapping_df.astype(str)
    return gene_mapping_df


def limit_rows_per_gene(df, stocksPerGene):
    return (
        df.groupby("ASSOCIATED_GENE")
        .apply(lambda x: x.sample(n=min(len(x), stocksPerGene)))
        .reset_index(drop=True)
    )


def affected_chromosomes(bloom_stock_no, bdsc_stocks):
    result = bdsc_stocks[bdsc_stocks["Stk #"] == bloom_stock_no][
        "All Affected Chromosomes"
    ]
    return result.iloc[0] if not result.empty else ""


def extract_affected_genes(genotype):
    if not isinstance(genotype, str):
        return "-"
    # Split genotype by delimiters and check each token against the set of gene symbols
    affected_genes = [gene for gene in all_gene_symbols if str(gene) in genotype]
    # Join the matches into a comma-separated string
    return ", ".join(map(str, affected_genes)) if affected_genes else "-"


def load_split_config(config_path):
    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def apply_filters(df, filters_config):
    filters = {}
    for filter_name, condition in filters_config.items():
        if condition["type"] == "contains":
            filters[filter_name] = df[condition["column"]].str.contains(
                condition["value"], case=True, na=False
            )
        elif condition["type"] == "doesnt_contain":
            filters[filter_name] = ~df[condition["column"]].str.contains(
                condition["value"], case=True, na=False
            )
        elif condition["type"] == "equals":
            filters[filter_name] = df[condition["column"]] == condition["value"]
        elif condition["type"] == "doesnt_equal":
            filters[filter_name] = df[condition["column"]] != condition["value"]
        # Add more conditions as needed
    return filters


def find_overlaps(stock_refs, gene_refs):
    if stock_refs == "-" or gene_refs == "-":
        return "-"
    
    # Split and strip references
    stock_set = set(ref.strip() for ref in stock_refs.split(", "))
    gene_set = set(ref.strip() for ref in gene_refs.split(", "))
    
    # Find overlap
    overlap = stock_set.intersection(gene_set)
    return ", ".join(sorted(overlap)) if overlap else "-"

# Function to extract overlapping terms and create 'GO_plus' and 'Relevant_GO_Terms'
def find_relevant_go_terms(row, relevant_terms):
    overlapping_terms = []

    # Convert relevant terms to lowercase for case-insensitive comparison
    relevant_terms_lower = {term.lower() for term in relevant_terms}

    # Process and clean each term in the row
    for term in row.split(','):
        term = term.strip()  # Remove whitespace
        if '~' in term:
            term = term.split('~')[-1]  # Remove everything before and including '~'
        elif ':' in term:
            term = term.split(':')[-1]  # Remove everything before and including ':'
        term_lower = term.lower()  # Convert term to lowercase

        # Check for an exact match in the relevant terms
        if term_lower in relevant_terms_lower:
            overlapping_terms.append(term_lower)

    # Return True/False and the unique overlapping relevant terms
    return pd.Series([bool(overlapping_terms), ', '.join(sorted(set(overlapping_terms)))])

def calculate_relevance_score(stock_sleep_papers, stock_paper_refs):
    if (
        stock_sleep_papers != "-"
    ):  # Allele sleep papers exist
        return 2
    elif (
        stock_paper_refs != "-"
    ):  # No overlaps, but STOCK_PAPER_REFS is not empty
        return 1
    else:  # STOCK_PAPER_REFS is empty
        return 0