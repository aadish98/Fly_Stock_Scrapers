# Script adds:
# - PMID
# - title/ abstract
# - and keyword analysis to csv files with FBrf's

# CMD Line args:
# - Directory Path,
# - Reference ID col,
# - Keywords csv filepath,
# - Keywords aggregate col

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd
import csv
import re

try:
    from Bio import Entrez
except ImportError:
    raise ImportError(
        "This script requires Biopython. Install via:\n"
        "    pip install biopython"
    )

# Set your email for Entrez queries
Entrez.email = "aadish98@gmail.com"

########################################
# A) Local Cache Configuration
########################################

# Folder and file path for your PMID -> title/abstract cache
PUBMED_CACHE_DIR = "/Users/AadishShah/Desktop/Projects.nosync/Allada Lab Downloads/Fly_Stock_Scrapers/Data/PubMed"
PUBMED_CACHE_FILE = os.path.join(PUBMED_CACHE_DIR, "pmid_to_title_abstract.csv")

# Global dictionary: pmid_cache["12345"] = {"title": "...", "abstract": "..."}
pmid_cache = {}

def ensure_cache_path_exists():
    """
    Create the PubMed cache directory if it doesn't already exist.
    """
    if not os.path.isdir(PUBMED_CACHE_DIR):
        os.makedirs(PUBMED_CACHE_DIR, exist_ok=True)

def aggregate_keywords(keywords):
    all_keywords = set()
    for entry in keywords:
        all_keywords.update(entry.split(", "))  # Split by ', ' and add to set (deduplication)

    all_keywords.discard("-")
    all_keywords = set(all_keywords)
    return ", ".join(sorted(all_keywords))  # Sort and join back to string

def load_pmid_cache():
    """
    Loads pmid_cache from PUBMED_CACHE_FILE if it exists.
    CSV columns: pmid, title, abstract
    """
    ensure_cache_path_exists()
    if not os.path.isfile(PUBMED_CACHE_FILE):
        print(f"[INFO] No existing cache found at {PUBMED_CACHE_FILE}. Starting fresh.")
        return

    cache_df = pd.read_csv(PUBMED_CACHE_FILE, dtype=str, keep_default_na=False)
    for _, row in cache_df.iterrows():
        pmid = row["pmid"]
        pmid_cache[pmid] = {
            "title": row["title"],
            "abstract": row["abstract"]
        }

    print(f"[INFO] Loaded {len(pmid_cache)} entries from local cache.")

def save_pmid_cache():
    """
    Saves the current pmid_cache to PUBMED_CACHE_FILE (CSV) with columns: pmid, title, abstract.
    """
    ensure_cache_path_exists()
    items = []
    for pmid, data in pmid_cache.items():
        items.append({
            "pmid": pmid,
            "title": data["title"],
            "abstract": data["abstract"]
        })
    df_cache = pd.DataFrame(items, columns=["pmid", "title", "abstract"])
    df_cache.to_csv(PUBMED_CACHE_FILE, index=False)
    print(f"[INFO] Saved updated cache with {len(df_cache)} PMIDs to {PUBMED_CACHE_FILE}.")

########################################
# B) FlyBase references
########################################

def load_flybase_refs(tsv_path):
    """
    Reads a FlyBase references TSV with columns:
      FBrf, PMID, PMCID, DOI, pub_type, miniref, pmid_added
    Returns a dict: { 'FBrf12345': ['1234567', ...], ... }
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df.fillna("", inplace=True)

    # Just keep FBrf and PMID
    df = df[["FBrf", "PMID"]].drop_duplicates()

    # Build a dictionary for easy FBrf -> PMIDs
    fbrf_to_pmid = {}
    for _, row in df.iterrows():
        fbrf = row["FBrf"].strip()
        pmid = row["PMID"].strip()
        if fbrf not in fbrf_to_pmid:
            fbrf_to_pmid[fbrf] = []
        if pmid:
            fbrf_to_pmid[fbrf].append(pmid)
    return fbrf_to_pmid

########################################
# C) PubMed fetch with caching
########################################

def get_or_fetch_pubmed_title_abstract(pmid):
    """
    Unified function that checks pmid_cache first:
      - If pmid is cached, return (title, abstract).
      - Otherwise, fetch from PubMed via esummary (for title) and efetch (for abstract),
        store in cache, then return it.
    """
    pmid = pmid.strip()
    if not pmid or pmid == "-":
        return "-", "-"

    # Check cache
    if pmid in pmid_cache:
        data = pmid_cache[pmid]
        return data["title"], data["abstract"]

    # Not cached; fetch fresh data.
    fetched_title = "-"
    fetched_abstract = "-"

    # 1) Fetch title via esummary
    try:
        summary_handle = Entrez.esummary(db="pubmed", id=pmid)
        summary_data = Entrez.read(summary_handle)
        summary_handle.close()

        if summary_data and isinstance(summary_data, list):
            doc_sum = summary_data[0]
            raw_title = doc_sum.get("Title", "")
            # Replace newlines with a placeholder if desired
            fetched_title = raw_title.replace("\n", " ").strip() if raw_title else "-"
    except Exception as e:
        print(f"[WARNING] Could not retrieve title for PMID={pmid}. Error: {e}")

    # 2) Fetch abstract via efetch
    try:
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        text = fetch_handle.read().strip()
        fetch_handle.close()
        if text:
            fetched_abstract = text.replace("\n", " ").strip()
    except Exception as e:
        print(f"[WARNING] Could not retrieve abstract for PMID={pmid}. Error: {e}")

    # Store in cache
    pmid_cache[pmid] = {
        "title": fetched_title,
        "abstract": fetched_abstract
    }

    return fetched_title, fetched_abstract

########################################
# D) Keyword matching
########################################

def load_keywords(keywords_csv):
    """
    Reads the keywords CSV, which must have a single column named 'keywords'.
    Returns a set of cleaned, non-empty keywords.
    """
    kw_df = pd.read_csv(keywords_csv, dtype=str).fillna("")
    if "keywords" not in kw_df.columns:
        raise ValueError("The keywords CSV must have a column named 'keywords'.")

    keywords_list = [k.strip() for k in kw_df["keywords"] if k.strip()]
    return set(keywords_list)

def find_matched_keywords(title_str, abstract_str, keywords):
    combined = (title_str + " " + abstract_str).lower()
    found = []

    for kw in keywords:
        kw_lower = kw.lower()
        # Match full words OR part of words at word boundaries
        pattern = rf"\b{re.escape(kw_lower)}\b"
        if re.search(pattern, combined):
            found.append(kw)

    return ", ".join(sorted(set(found))) if found else "-"


########################################
# E) Main logic (same overall flow)
########################################

def main():
    parser = argparse.ArgumentParser(
        description="""
        1) Reads a FlyBase references TSV (hard-coded path),
        2) For each CSV in a directory (skipping 'summary' in filename),
           uses a user-specified column containing one or more FBrfs (split by
           commas or newlines) to look up PMIDs, fetch PubMed titles/abstracts
           (with local caching), and match keywords from a CSV.
        3) Overwrites each CSV, adding columns 'PMID', 'title', 'abstract', and the
           user-specified 'colname' for matched keywords.
        """
    )
    parser.add_argument("directory", help="Path to folder containing CSV files.")
    parser.add_argument("ref_col", help="Name of the column in the CSV that has FBrfs (may contain multiple, separated by commas/newlines).")
    parser.add_argument("keywords_csv", help="Path to CSV with a single 'keywords' column.")
    parser.add_argument("colname", help="Name of the new column storing matched keywords (e.g. 'keywords_in_title_abstract').")
    args = parser.parse_args()

    # 1) Load local PubMed cache (or start fresh if none)
    load_pmid_cache()

    # 2) Hard-coded path to FlyBase references TSV
    flybase_refs_path = (
        "/Users/AadishShah/Desktop/Projects.nosync/"
        "Allada Lab Downloads/Fly_Stock_Scrapers/Data/FlyBase/FlyBase_References/"
        "fbrf_pmid_pmcid_doi_fb_2024_05.tsv"
    )
    fbrf_dict = load_flybase_refs(flybase_refs_path)

    # 3) Load keywords
    keywords = load_keywords(args.keywords_csv)

    # 4) Process each CSV in the directory
    csv_files = glob.glob(os.path.join(args.directory, "*.csv"))
    if not csv_files:
        print(f"[INFO] No CSV files found in {args.directory}.")
        return

    for csv_path in csv_files:
        basename = os.path.basename(csv_path)
        if "summary" in basename or "keyword" in basename:
            print(f"[INFO] Skipping file with 'summary' or 'keyword' in name: {csv_path}")
            continue

        print(f"\n[INFO] Processing: {csv_path} ...")
        df = pd.read_csv(csv_path, dtype=str).fillna("")

        # Check if user-specified column exists
        if args.ref_col not in df.columns:
            print(f"[WARNING] Column '{args.ref_col}' not found in {csv_path}. Skipping.")
            continue

        df[args.ref_col] = df[args.ref_col].str.split(", ")  # Split by ", " into lists

        # Now, explode the column
        df["ref_id_placeholder"] = df[args.ref_col]
        df = df.explode("ref_id_placeholder")
        # We'll create new lists to store final strings for PMID, title, abstract, and keywords
        pmid_col_vals = []
        title_col_vals = []
        abstract_col_vals = []
        keyword_col_vals = []

        # For each row, parse the references from df[args.ref_col].
        # This column may have multiple references separated by commas or newlines.
        for _, row in df.iterrows():
            ref_string = row['ref_id_placeholder'].strip()
            if not ref_string:
                # No FBrfs
                pmid_col_vals.append("-")
                title_col_vals.append("-")
                abstract_col_vals.append("-")
                keyword_col_vals.append("-")
                continue

            # Split by commas/newlines/semicolons
            refs_list = [r.strip() for r in re.split(r'[,\n;]+', ref_string) if r.strip()]

            # Now gather all PMIDs for all FBrfs
            row_pmids = []
            for fbrf in refs_list:
                if fbrf in fbrf_dict:
                    row_pmids.extend(fbrf_dict[fbrf])
                else:
                    row_pmids.append("-")

            # De-duplicate pmids, if desired
            row_pmids = list({p for p in row_pmids if p})
            if not row_pmids:
                row_pmids = ["-"]

            # Join them with semicolon for storage
            pmid_str = "; ".join(row_pmids)
            pmid_col_vals.append(pmid_str)

            # For each PMID, fetch from cache or PubMed, then join with semicolon
            row_titles = []
            row_abstracts = []
            for pmid_val in row_pmids:
                if pmid_val == "-":
                    row_titles.append("-")
                    row_abstracts.append("-")
                else:
                    t, a = get_or_fetch_pubmed_title_abstract(pmid_val)
                    row_titles.append(t)
                    row_abstracts.append(a)

            title_str = "; ".join(row_titles)
            abstract_str = "; ".join(row_abstracts)

            title_col_vals.append(title_str)
            abstract_col_vals.append(abstract_str)

            # Match keywords
            matched = find_matched_keywords(title_str, abstract_str, keywords)
            keyword_col_vals.append(matched)

        # Insert new columns
        df["PMID"] = pmid_col_vals
        df["title"] = title_col_vals
        df["abstract"] = abstract_col_vals
        df[args.colname] = keyword_col_vals

        # Reorder columns
        original_cols = [c for c in df.columns if c not in ("PMID", "title", "abstract", args.colname)]
        final_cols = original_cols + ["PMID", "title", "abstract", args.colname, "all_aggregated_gene_ref_keywords", "total_count_of_gene_keywords"]
        aggregated_keywords = df.groupby("ext_gene")[args.colname].agg(aggregate_keywords)
        keyword_counts = aggregated_keywords.str.split(", ").apply(len)
        # Merge aggregated column back to original DataFrame
        df["all_aggregated_gene_ref_keywords"] = df["ext_gene"].map(aggregated_keywords)
        df["total_count_of_gene_keywords"] = df["ext_gene"].map(keyword_counts)
        df = df[final_cols]

        # Overwrite CSV with quoting
        df.to_csv(csv_path, index=False, quoting=csv.QUOTE_ALL)
        print(f"[INFO] Updated CSV saved: {csv_path}")

    # 5) Save updated cache
    save_pmid_cache()

if __name__ == "__main__":
    main()
