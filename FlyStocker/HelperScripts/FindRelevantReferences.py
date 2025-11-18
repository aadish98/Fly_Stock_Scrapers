#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd
import re
from tqdm import tqdm  #  Import tqdm for progress bars
import csv

def contains_keyword(text, keyword):
    """
    Checks whether 'keyword' appears in 'text' as a distinct substring,
    i.e. not as part of a longer word.
    """
    # text = text.lower()
    # keyword = keyword.lower()
    # start = 0
    # while True:
    #     pos = text.find(keyword, start)
    #     if pos == -1:
    #         return False
    #     # Check the character before and after the keyword (if they exist)
    #     before_ok = (pos == 0) or (not text[pos - 1].isalpha())
    #     after_index = pos + len(keyword)
    #     after_ok = (after_index == len(text)) or (not text[after_index].isalpha())
    #     if before_ok and after_ok:
    #         return True
    #     # Otherwise, look for another occurrence
    #     start = pos + 1
    return keyword.lower() in text.lower()

def find_matched_keywords(title, abstract, keywords_set):
    """
    Checks if any keyword from keywords_set appears in the given title + abstract
    as a distinct word (i.e. not as part of another word). Returns a unique,
    alphabetically sorted list of matched keywords.
    """
    combined_text = f"{title} {abstract}"
    matched_keywords = [kw for kw in keywords_set if contains_keyword(combined_text, kw)]
    return "; ".join(sorted(set(matched_keywords))) if matched_keywords else None


def aggregate_keywords(df):
    """
    Aggregates 'keywords_in_title_abstract' per unique 'flybase_gene_id',
    ensuring the keywords are unique and alphabetically sorted.
    """

    def sort_keywords(keyword_string):
        # Split by semicolon, remove extra spaces, filter out empty values,
        # then deduplicate and sort.
        keywords = [kw.strip() for kw in keyword_string.split(";") if kw.strip()]
        return "; ".join(sorted(set(keywords)))

    aggregated_keywords = (
        df.groupby("flybase_gene_id")["keywords_in_title_abstract"]
        .apply(lambda x: sort_keywords("; ".join(x.dropna())))
    )

    return df.assign(all_aggregated_genes_refs_keywords=df["flybase_gene_id"].map(aggregated_keywords))


def process_files(refs_dir, keywords_csv, save_dir):
    """
    Processes each CSV in the REFS directory by searching for keywords
    in the Title and Abstract, adding the two new columns, and saving the updated files
    in a subfolder inside REFS, with a progress bar.
    Also generates a summary of keyword occurrences across all files, split by individual keywords.
    """
    # Load keyword list
    keywords = pd.read_csv(keywords_csv, dtype=str).fillna("")
    if "keywords" not in keywords.columns:
        raise ValueError("The keywords CSV must have a column named 'keywords'.")

    keywords_set = set(keywords["keywords"].str.strip())

    # Create output directory
    output_dir = os.path.join(refs_dir, save_dir)
    os.makedirs(output_dir, exist_ok=True)
    print(f"[INFO] Saving processed files in: {output_dir}")

    # Find all CSV files in REFS directory
    csv_files = glob.glob(os.path.join(refs_dir, "*.csv"))
    if not csv_files:
        print(f"[ERROR] No CSV files found in {refs_dir}. Exiting.")
        return

    #  Progress bar for files
    summary_list = []

    for csv_file in tqdm(csv_files, desc="Processing CSV files", unit="file"):
        summary_dict = {}
        df = pd.read_csv(csv_file, dtype=str, keep_default_na=False,
                         quoting=csv.QUOTE_ALL,  # Ensure proper handling of quoted fields
                         lineterminator="\n")  # Prevent unintended line breaks
        summary_dict["Filename"] = os.path.basename(csv_file)

        # Check if required columns exist
        if not {"flybase_gene_id", "title", "abstract", "flybase_references"}.issubset(df.columns):
            print(f"[WARNING] Required columns missing in {csv_file}. Skipping.")
            continue

        tqdm.pandas(desc="Matching Keywords")
        df["keywords_in_title_abstract"] = df.progress_apply(
            lambda row: find_matched_keywords(row["title"], row["abstract"], keywords_set), axis=1
        )

        df.dropna(subset=["keywords_in_title_abstract"], inplace=True)

        df = aggregate_keywords(df)

        #  Initialize dictionary to track keyword-specific counts
        keyword_gene_counts = {kw: set() for kw in keywords_set}  # Store unique genes per keyword
        keyword_ref_counts = {kw: set() for kw in keywords_set}  # Store unique references per keyword

        #  Populate keyword-specific counts
        for _, row in df.iterrows():
            matched_keywords = row["keywords_in_title_abstract"].split("; ")
            gene_id = row["flybase_gene_id"]
            ref_id = row["flybase_references"]

            for kw in matched_keywords:
                if kw in keyword_gene_counts:
                    keyword_gene_counts[kw].add(gene_id)  # Track unique genes per keyword
                if kw in keyword_ref_counts:
                    keyword_ref_counts[kw].add(ref_id)  # Track unique references per keyword

        #  Store the number of unique hits per keyword
        summary_dict["Genes with keyword hits: "] = df.flybase_gene_id.nunique()

        summary_dict["Refs with keyword hits: "] = df.flybase_references.nunique()
        for kw in keywords_set:
            summary_dict[f"Genes with {kw}"] = len(keyword_gene_counts[kw])
            summary_dict[f"Refs with {kw}"] = len(keyword_ref_counts[kw])

        #  Append to summary
        summary_list.append(summary_dict)

        #  Save updated file
        output_file = os.path.join(output_dir, os.path.basename(csv_file))
        if "frequency" in df.columns:
            df = df.sort_values(by=["frequency", "flybase_gene_id"], ascending=[False, True])
        df.to_csv(output_file, index=False, quoting=csv.QUOTE_ALL)
        print(f"[INFO] Updated file saved: {output_file}")

    summary_df = pd.DataFrame(summary_list)


    #  Save the summary as CSV
    summary_df.to_csv(os.path.join(output_dir, "summary_df.csv"), index=False, encoding="utf-8", quoting=csv.QUOTE_ALL)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Match keywords in Title & Abstract and aggregate by flybase_gene_id.")
    parser.add_argument("refs_dir", help="Path to REFS directory containing processed CSV files.")
    parser.add_argument("keywords_csv", help="Path to CSV containing keywords (must have 'keywords' column).")
    parser.add_argument("save_dir", help="Subdirectory inside REFS to save processed files.")
    args = parser.parse_args()

    process_files(args.refs_dir, args.keywords_csv, args.save_dir)
