#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import time
import csv

# Set your email for NCBI Entrez API
Entrez.email = "aadish98@gmail.com"
pd.set_option("display.max_colwidth", None)  #  Prevents truncation in pandas display

# --- PubMed Cache Configuration ---
PUBMED_CACHE_DIR = "../../Data/PubMed"
PUBMED_CACHE_FILE = os.path.join(PUBMED_CACHE_DIR, "pmid_to_title_abstract.csv")


def load_pmid_cache():
    """Load cached PMIDs from PUBMED_CACHE_FILE."""
    if not os.path.exists(PUBMED_CACHE_FILE):
        print(f"[INFO] No existing PubMed cache found. Starting fresh.")
        return {}

    cache_df = pd.read_csv(PUBMED_CACHE_FILE, dtype=str, keep_default_na=False, quoting=csv.QUOTE_ALL)
    cached_data = {row["pmid"]: {"title": row["title"], "abstract": row["abstract"]} for _, row in cache_df.iterrows()}

    print(f"[INFO] Loaded {len(cached_data)} entries from PubMed cache.")
    return cached_data


def save_pmid_cache(new_data):
    """Append new PMIDs to PUBMED_CACHE_FILE."""
    if not new_data:
        return

    new_df = pd.DataFrame(new_data).T.reset_index().rename(columns={"index": "pmid"})
    if os.path.exists(PUBMED_CACHE_FILE):
        existing_df = pd.read_csv(PUBMED_CACHE_FILE, dtype=str, keep_default_na=False, quoting=csv.QUOTE_ALL)
        updated_df = pd.concat([existing_df, new_df], ignore_index=True).drop_duplicates(subset=["pmid"])
    else:
        updated_df = new_df

    os.makedirs(PUBMED_CACHE_DIR, exist_ok=True)
    updated_df.to_csv(PUBMED_CACHE_FILE, index=False)
    print(f"[INFO] Saved {len(new_data)} new PMIDs to cache.")


def fetch_pubmed_data(pmids, pmid_cache):
    """Fetch Title & Abstract for PMIDs not in cache using NCBI Entrez API."""

    batch_size = 200
    new_data = {}

    pmids_to_fetch = [pmid for pmid in pmids if pmid not in pmid_cache]

    if not pmids_to_fetch:
        print("[INFO] All PMIDs found in cache. No need to query PubMed.")
        return pmid_cache, new_data  #  No new PMIDs were fetched

    for i in tqdm(range(0, len(pmids_to_fetch), batch_size), desc="Fetching PubMed data"):
        batch_pmids = pmids_to_fetch[i: i + batch_size]  #  Correctly slicing batch
        pmid_str = ",".join(batch_pmids)  #  Join PMIDs properly

        try:
            #  Fetch article metadata in MEDLINE format
            handle = Entrez.efetch(db="pubmed", id=pmid_str, rettype="medline", retmode="text")
            records = handle.read().strip().split("\n\n")  #  Ensure no extra newlines
            handle.close()

            for record in records:
                lines = record.strip().split("\n")

                pmid, title, abstract = None, None, None  #  Keep None for missing fields
                current_field = None  #  Track which field we're processing

                for line in lines:
                    line = line.strip()

                    if line.startswith("PMID- "):
                        pmid = line.replace("PMID- ", "").strip()

                    elif line.startswith("TI  - "):  #  Start of title
                        title = line.replace("TI  - ", "").strip()
                        current_field = "title"

                    elif line.startswith("AB  - "):  #  Start of abstract
                        abstract = line.replace("AB  - ", "").strip()
                        current_field = "abstract"

                    elif line.startswith("OT  - "):  #  Stop appending if other fields start
                        current_field = None

                    elif current_field == "title":
                        title += " " + line.strip()  #  Append continuation lines to title

                    elif current_field == "abstract":
                        abstract += " " + line.strip()  #  Append continuation lines to abstract

                #  Ensure pmid is valid before storing
                if pmid and pmid.isdigit():
                    new_data[pmid] = {"title": title.strip() if title else None,
                                      "abstract": abstract.strip() if abstract else None}

        except Exception as e:
            print(f"[ERROR] Failed to fetch batch: {batch_pmids}. Error: {e}")
            time.sleep(10)  #  Increased sleep time to handle rate limits
            continue  #  Continue processing the next batch even if one fails

    #  Merge new data into cache
    pmid_cache.update(new_data)
    return pmid_cache, new_data

def main():
    parser = argparse.ArgumentParser(
        description="Left-join each CSV in a folder with a filtered FlyBase references TSV.")
    parser.add_argument("directory", help="Path to the directory containing CSV files.")
    args = parser.parse_args()

    input_dir = args.directory
    output_dir = os.path.join(input_dir, "REFS")  # Create output directory path

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # --- Load PubMed Cache ---
    pmid_cache = load_pmid_cache()

    # --- Load FlyBase references ---
    refs_path = (
        "../../Data/FlyBase/FlyBase_References/"
        "entity_publication_fb_2024_06.tsv"
    )

    refs_type_path = (
        "../../Data/FlyBase/FlyBase_References/"
        "fbrf_pmid_pmcid_doi_fb_2024_06.tsv"
    )

    blacklist_path = (
        "../../Data/FlyBase/custom/blacklisted_refs.csv"
    )

    # --- Load FlyBase gene-to-reference mapping ---
    refs_gene_df = pd.read_csv(refs_path, sep='\t', dtype=str, low_memory=False, quoting=csv.QUOTE_ALL)
    refs_gene_df = refs_gene_df[refs_gene_df["entity_id"].str.startswith("FBgn", na=False)]
    refs_gene_df = refs_gene_df[["FlyBase_publication_id", "entity_id"]]
    refs_gene_df.rename(columns={"FlyBase_publication_id": "FBrf"}, inplace=True)

    # --- Load reference types and filter for 'paper' only ---
    refs_type_df = pd.read_csv(refs_type_path, sep='\t', dtype=str, low_memory=False, quoting=csv.QUOTE_ALL)
    refs_type_df = refs_type_df[refs_type_df["pub_type"] == "paper"]
    refs_type_df = refs_type_df[["FBrf", "PMID", "PMCID", "DOI"]]


    # --- Load blacklisted references ---
    blacklist_df = pd.read_csv(blacklist_path, sep=',', dtype=str, low_memory=False)
    blacklisted_refs = set(blacklist_df["FBrf"].dropna())

    # Remove blacklisted references from refs_gene_df
    refs_gene_df = refs_gene_df[~refs_gene_df["FBrf"].isin(blacklisted_refs)]

    # Merge refs_gene_df with refs_type_df to keep only references of type 'paper' and retain PMIDs
    refs_gene_df = refs_gene_df.merge(refs_type_df, how="inner", on="FBrf")

    # --- Aggregate references and PMIDs per gene ---
    gene_reference_map = refs_gene_df.groupby("entity_id")["FBrf"].agg("; ".join).to_dict()
    gene_pmid_map = refs_gene_df.groupby("entity_id")["PMID"].agg("; ".join).to_dict()

    # --- Fetch Titles & Abstracts from PubMed ---
    all_pmids = set(refs_gene_df["PMID"].dropna())
    pmid_cache, new_pmid_data = fetch_pubmed_data(all_pmids, pmid_cache)

    # Save new PMIDs to cache
    save_pmid_cache(new_pmid_data)

    # --- Process each CSV file in the input directory ---
    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    if not csv_files:
        print(f"[INFO] No CSV files found in {input_dir}. Exiting.")
        return

    for csv_file in csv_files:
        print(f"\n[INFO] Processing {csv_file} ...")
        df = pd.read_csv(csv_file, dtype=str, keep_default_na=False, quoting=csv.QUOTE_ALL, lineterminator="\n")

        if "flybase_gene_id" not in df.columns:
            print(f"[WARNING] 'flybase_gene_id' column not found in {csv_file}, skipping.")
            continue

        df["flybase_references"] = df["flybase_gene_id"].map(gene_reference_map).fillna("")
        df["PMID"] = df["flybase_gene_id"].map(gene_pmid_map).fillna("")

        df["flybase_references"] = df["flybase_references"].str.split("; ")
        df["PMID"] = df["PMID"].str.split("; ")
        df = df.explode(["flybase_references", "PMID"]).reset_index(drop=True)

        df["title"] = df["PMID"].map(lambda pmid: pmid_cache.get(pmid, {}).get("title", ""))
        df["abstract"] = df["PMID"].map(lambda pmid: pmid_cache.get(pmid, {}).get("abstract", ""))

        df = df.merge(refs_gene_df[["FBrf",	"PMCID",	"DOI"]], how = "left", left_on = "flybase_references", right_on = "FBrf")
        df = df.dropna(subset = ["flybase_references"])
        df = df.drop(columns = ["FBrf"])
        df = df[["flybase_references","flybase_gene_id", "ext_gene", "PMID","PMCID", "DOI", "title", "abstract"]]
        df = df[df.flybase_references != ""]
        df["title_populated"] = df["title"].apply(lambda x: 1 if isinstance(x, str) and x.strip() else 0)
        df["abstract_populated"] = df["abstract"].apply(lambda x: 1 if isinstance(x, str) and x.strip() else 0)

        # Create a priority column: prioritize rows with both title and abstract populated, followed by title-only, abstract-only, then none
        df["priority"] = df["title_populated"] + df["abstract_populated"]

        # Sort DataFrame so that rows with populated title/abstract are prioritized
        df = df.sort_values(by=["flybase_references", "flybase_gene_id", "priority"], ascending=[True, True, False])

        # Drop duplicates, keeping the highest-priority row for each (flybase_references, flybase_gene_id) combination
        df = df.drop_duplicates(subset=["flybase_references", "flybase_gene_id"])

        # Drop the helper columns used for prioritization
        df = df.drop(columns=["title_populated", "abstract_populated", "priority"])
        output_file = os.path.join(output_dir, os.path.basename(csv_file))
        df.to_csv(output_file, index=False, encoding="utf-8", quoting=csv.QUOTE_ALL)

        missing_pmcid_doi_count = df[df["PMCID"].isna()].shape[0]

        print(f"[INFO] % of Rows with missing PMCID and DOI: {missing_pmcid_doi_count}/{df.shape[0]}")

        print(f"[INFO] Saved processed file to: {output_file}")


if __name__ == "__main__":
    main()
