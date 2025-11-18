# Script adds publications (ref_id = column of FBrf's):
#   -DOI,
#   -Title,
#   -miniref (author, journal, year),
#   -and whether its a 'blacklisted' ref (i.e. > 50 genes)

# Input
# - CMD LINE:
# - Directory path with csv files

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import glob
import argparse
from crossref_commons.retrieval import get_publication_as_json

########################################
# Configuration
########################################

# Path to your persistent CSV cache
CACHE_FILE = "/Users/AadishShah/Desktop/Projects.nosync/Allada Lab Downloads/Fly_Stock_Scrapers/Data/CrossRefs/doi_title_cache.csv"

# Path to your blacklisted references
BLACKLIST_FILE = "../../Data/FlyBase/custom/blacklisted_refs.csv"

# Path(s) to your FlyBase tsv data (FBrf -> DOI, etc.)
FLYBASE_FOLDER = "../../Data/FlyBase/"
FLYBASE_TSV_PATTERN = "fbrf_pmid_*.tsv"

# Columns we will create (you can adjust names to your liking)
DOI_COL            = "PAPER_DOI"
TITLE_COL          = "PAPER_TITLE"
MINIREF_COL        = "PAPER_MINIREF"
BLACKLISTED_COL    = "BLACKLISTED_REF"  # new column indicating blacklisted status

# A column in your CSV(s) that contains references (FBrf IDs).
# If your CSV references are stored differently, adapt parsing accordingly.
REFERENCE_COL = "ref_id"

########################################
# Step 1: Load or initialize the cache
########################################

doi_metadata_cache = {}

if os.path.exists(CACHE_FILE):
    cache_df = pd.read_csv(CACHE_FILE, keep_default_na=False)
    for _, row in cache_df.iterrows():
        doi = row["doi"]
        doi_metadata_cache[doi] = {
            "title": row["title"]
        }
    print(f"[INFO] Loaded {len(doi_metadata_cache)} entries from cache.")
else:
    print(f"[INFO] No existing cache found at {CACHE_FILE}; starting fresh.")

########################################
# Step 2: Function to fetch Crossref metadata
########################################

def fetch_metadata_from_crossref(doi):
    """
    Retrieves title from Crossref.
    Uses a global cache to avoid repeated calls.
    """
    if not doi or doi.strip() == "-":
        return {"title": "-"}

    doi = doi.strip()
    if doi in doi_metadata_cache:
        return doi_metadata_cache[doi]

    try:
        paper_json = get_publication_as_json(doi)
        title = paper_json.get("title", ["Title not found"])[0]
        result = {"title": title}
    except Exception as e:
        result = {"title": f"Error fetching: {e}"}

    # Store in cache
    doi_metadata_cache[doi] = result
    return result

########################################
# Step 3: Load FlyBase paper data
########################################

def load_flybase_paper_data():
    tsv_files = glob.glob(os.path.join(FLYBASE_FOLDER, "**", FLYBASE_TSV_PATTERN), recursive=True)
    if not tsv_files:
        raise FileNotFoundError(f"No TSV files found in {FLYBASE_FOLDER} pattern '{FLYBASE_TSV_PATTERN}'")

    dfs = []
    for file_path in tsv_files:
        df = pd.read_csv(file_path, delimiter='\t', low_memory=False)
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["FBrf"])
    if "pub_type" in combined_df.columns:
        combined_df = combined_df[combined_df["pub_type"] == "paper"]

    # We'll keep only FBrf, DOI, and miniref
    keep_cols = ["FBrf", "DOI", "miniref"]
    keep_cols = [c for c in keep_cols if c in combined_df.columns]
    combined_df = combined_df[keep_cols].drop_duplicates()

    if "FBrf" not in combined_df.columns:
        raise ValueError("Combined FlyBase data has no 'FBrf' column, cannot proceed.")

    combined_df.set_index("FBrf", inplace=True)
    return combined_df

flybase_paper_df = load_flybase_paper_data()

########################################
# Step 4: Load blacklisted references
########################################

if os.path.exists(BLACKLIST_FILE):
    blacklisted_refs = pd.read_csv(BLACKLIST_FILE)
    blacklisted_refs_set = set(blacklisted_refs["FBrf"])
else:
    print(f"[WARNING] Blacklist file not found at {BLACKLIST_FILE}. No references will be blacklisted.")
    blacklisted_refs_set = set()

########################################
# Step 5: Functions to parse references, convert to DOIs & minirefs
########################################

def get_doi_for_fbrf(fbrf_id):
    if fbrf_id in flybase_paper_df.index:
        doi_value = flybase_paper_df.loc[fbrf_id, "DOI"]
        if isinstance(doi_value, pd.Series):
            doi_value = doi_value.dropna().iloc[0] if not doi_value.dropna().empty else "-"
        return doi_value if pd.notna(doi_value) else "-"
    return "-"

def get_miniref_for_fbrf(fbrf_id):
    if fbrf_id in flybase_paper_df.index and "miniref" in flybase_paper_df.columns:
        miniref_value = flybase_paper_df.loc[fbrf_id, "miniref"]
        if isinstance(miniref_value, pd.Series):
            miniref_value = miniref_value.dropna().iloc[0] if not miniref_value.dropna().empty else "-"
        return miniref_value if pd.notna(miniref_value) else "-"
    return "-"

def parse_and_fetch_fbrf_refs(ref_string):
    """
    Returns three parallel lists (same length):
    - A list of DOIs
    - A list of minirefs
    - A list of blacklist indicators ("YES"/"NO")
    """
    if not isinstance(ref_string, str) or not ref_string.strip() or ref_string.strip() == "-":
        return ["-"], ["-"], ["-"]

    fbrf_list = [r.strip() for r in ref_string.split(",") if r.strip()]
    doi_list = []
    miniref_list = []
    blacklisted_list = []

    for fbrf_id in fbrf_list:
        doi_list.append(get_doi_for_fbrf(fbrf_id))
        miniref_list.append(get_miniref_for_fbrf(fbrf_id))

        # Mark as "YES" if it's in the blacklisted set, else "NO"
        if fbrf_id in blacklisted_refs_set:
            blacklisted_list.append("YES")
        else:
            blacklisted_list.append("NO")

    if not doi_list:
        return ["-"], ["-"], ["-"]
    return doi_list, miniref_list, blacklisted_list

########################################
# Step 6: Main logic
########################################

def main():
    parser = argparse.ArgumentParser(description="Add Crossref metadata (DOI, title, miniref, blacklisted) columns to CSV files in a folder.")
    parser.add_argument("folder_path", help="Path to the folder containing CSV files.")
    args = parser.parse_args()

    input_folder = args.folder_path
    if not os.path.isdir(input_folder):
        print(f"[ERROR] The path '{input_folder}' is not a valid directory.")
        return

    csv_files = glob.glob(os.path.join(input_folder, "*.csv"))

    if not csv_files:
        print(f"[INFO] No CSV files found in {input_folder}. Exiting.")
        return

    for csv_file in csv_files:
        print(f"\n[INFO] Processing {csv_file} ...")
        df = pd.read_csv(csv_file, keep_default_na=False)

        if REFERENCE_COL not in df.columns:
            print(f"[WARNING] Column '{REFERENCE_COL}' not found in {csv_file}. Skipping file.")
            continue

        # Build new columns
        all_dois = []
        all_titles = []
        all_minirefs = []
        all_blacklisted = []

        for ref_str in df[REFERENCE_COL]:
            # 1) Parse references into DOIs, minirefs, and blacklisted flags
            doi_list, miniref_list, blacklisted_list = parse_and_fetch_fbrf_refs(ref_str)

            # 2) For each DOI, fetch metadata (title) from Crossref
            titles = []
            for doi in doi_list:
                meta = fetch_metadata_from_crossref(doi)
                titles.append(meta["title"])

            # Combine multiple references' metadata with line breaks
            all_dois.append("\n".join(doi_list))
            all_titles.append("\n".join(titles))
            all_minirefs.append("\n".join(miniref_list))
            all_blacklisted.append("\n".join(blacklisted_list))

        # Insert or update columns
        df[DOI_COL]            = all_dois
        df[TITLE_COL]          = all_titles
        df[MINIREF_COL]        = all_minirefs
        df[BLACKLISTED_COL]    = all_blacklisted

        df.to_csv(csv_file, index=False)
        print(f"[INFO] Updated file saved: {csv_file}")

    save_cache()

def save_cache():
    items = []
    for doi, meta in doi_metadata_cache.items():
        items.append({
            "doi":   doi,
            "title": meta["title"]
        })
    cache_df = pd.DataFrame(items, columns=["doi", "title"])
    cache_df.to_csv(CACHE_FILE, index=False)
    print(f"[INFO] Saved updated cache of {len(doi_metadata_cache)} DOIs to {CACHE_FILE}.")

if __name__ == "__main__":
    main()
