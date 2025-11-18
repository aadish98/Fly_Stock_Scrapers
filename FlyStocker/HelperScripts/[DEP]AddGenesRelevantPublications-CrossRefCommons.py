# Script adds genes:
# - Aggregate (relevant) publications (FBrf's),
# - Aggregate publication titles,
# - Aggregate pub. miniref
# Script also adds a summary table of % relevant publication coverage over each gene set in the directory

# CMD Line args: - Directory path,
# - (flags) sleep, neurodegen, neurodegen_specific, TDP43,
# - Sortby parameter


import pandas as pd
import sys
import os
import glob
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Core/')))
from AddGeneMetadata import (
    get_sleep_genes, get_neurodegen_genes, get_neurodegen_specific_genes, get_tdp43_genes,
    get_sleep_paper_refs, get_nd_paper_refs, get_nd_specific_paper_refs, get_tdp43_paper_refs
)
import argparse
import numpy as np
import subprocess
import crossref_commons.retrieval
from crossref_commons.retrieval import get_publication_as_json
from openpyxl.utils import get_column_letter
from openpyxl.styles import Alignment
import math
MAX_COL_WIDTH = 100
LINE_HEIGHT    = 15 # ~points per line, adjust as needed
# Path to your persistent CSV cache
CACHE_FILE = "/Users/AadishShah/Desktop/Projects.nosync/Allada Lab Downloads/Fly_Stock_Scrapers/Data/CrossRefs/doi_title_cache.csv"

########################################
# Step 1: Load the existing cache if it exists
########################################

doi_title_cache = {}  # global cache

if os.path.exists(CACHE_FILE):
    cache_df = pd.read_csv(CACHE_FILE, keep_default_na=False)
    # Convert the DataFrame to a dict mapping { doi: title }
    doi_title_cache = dict(zip(cache_df["doi"], cache_df["title"]))
    print(f"[INFO] Loaded {len(doi_title_cache)} entries from cache.")
else:
    doi_title_cache = {}
    print("[INFO] No existing cache found; starting fresh.")



def fetch_title_from_crossref(doi):
    if doi in doi_title_cache:
        return doi_title_cache[doi]
    try:
        paper = get_publication_as_json(doi)
        if "title" in paper and paper["title"]:
            result = paper["title"][0]
        else:
            result = "Title not found"
    except Exception as e:
        result = f"Error: {e}"
    doi_title_cache[doi] = result
    return result

def fetch_titles_for_dois(doi_string):
    if not isinstance(doi_string, str) or doi_string.strip() == '-' or not doi_string.strip():
        return '-'
    doi_list = [d.strip() for d in doi_string.split('\n') if d.strip()]
    title_list = [fetch_title_from_crossref(doi) for doi in doi_list]
    return "\n".join(title_list)

def get_ref_info(ref_string, fbrf_paper_df):
    if not isinstance(ref_string, str) or ref_string.strip() == '-':
        return "-", "-"
    ref_list = [r.strip() for r in ref_string.split(',') if r.strip()]
    doi_list = []
    miniref_list = []
    for ref in ref_list:
        if ref in fbrf_paper_df.index:
            row = fbrf_paper_df.loc[ref]
            doi_value = str(row["DOI"]) if pd.notna(row["DOI"]) else "-"
            miniref_value = str(row["miniref"]) if pd.notna(row["miniref"]) else "-"
            doi_list.append(doi_value)
            miniref_list.append(miniref_value)
        else:
            doi_list.append("-")
            miniref_list.append("-")
    return "\n".join(doi_list), "\n".join(miniref_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("target_dir")
    parser.add_argument("-neurodegen", action="store_true")
    parser.add_argument("-neurodegen_specific", action="store_true")
    parser.add_argument("-sleep", action="store_true")
    parser.add_argument("-TDP43", action="store_true")
    parser.add_argument("-sortby", choices=["frequency"])
    args = parser.parse_args()

    target_dir = args.target_dir
    neurodegen_flag = args.neurodegen
    sleep_flag = args.sleep
    sortby = args.sortby

    # 1) For efficiency, call GetFBgnIDs.py once (if needed)
    #    if you must do it per-file, keep your code as is.
    #    Or skip if you only need it once total.
    subprocess.run(["python3", "GetFBgnIDs.py", target_dir, "ext_gene"],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    # 2) Load blacklists & fbrf paper DataFrame once
    blacklisted_refs = pd.read_csv("../../Data/FlyBase/custom/blacklisted_refs.csv")
    blacklisted_refs_set = set(blacklisted_refs["FBrf"])

    folder_path = '../../Data/FlyBase/'
    tsv_files = glob.glob(os.path.join(folder_path, '**/fbrf_pmid_*.tsv'), recursive=True)
    flybase_dict = {}
    for file in tsv_files:
        file_name = os.path.splitext(os.path.basename(file))[0]
        df = pd.read_csv(file, delimiter='\t', low_memory=False)
        flybase_dict[file_name] = df
    fbrf_paper_df = flybase_dict["fbrf_pmid_pmcid_doi_fb_2024_05"]
    fbrf_paper_df = fbrf_paper_df.loc[fbrf_paper_df['pub_type'] == 'paper', ['FBrf', 'DOI', 'miniref']]
    fbrf_paper_df.set_index('FBrf', inplace=True)

    # 3) Process each CSV in the target_dir
    target_files = glob.glob(os.path.join(target_dir, "**.csv"), recursive=True)
    summary_df_temp = []
    for target_file in target_files:
        if "summary" in target_file or "single" in target_file:
            continue

        meta_df = pd.read_csv(target_file, keep_default_na=False)
        summary_dict = {"Filename": os.path.basename(target_file)}
        print(os.path.basename(target_file))
        print(f'Number of total genes: {meta_df.shape[0]}')

        summary_dict["Total Genes"] = meta_df["flybase_gene_id"].nunique() if "flybase_gene_id" in meta_df.columns else 0

        if sleep_flag and "flybase_gene_id" in meta_df.columns:
            meta_df["GENE_HAS_SLEEP_PAPER_REF"]  = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_sleep_genes(fbgnid)
            )
            meta_df["GENE_SLEEP_PAPER_REFS"]    = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_sleep_paper_refs(fbgnid)
            )
            pct_sleep = (meta_df[meta_df["GENE_HAS_SLEEP_PAPER_REF"]]["flybase_gene_id"].nunique() / meta_df["flybase_gene_id"].nunique() * 100)
            print(f"Percentage of genes with sleep reference: {pct_sleep:.2f}%")
            summary_dict["% of genes with sleep paper ref"] = f"{pct_sleep:.2f}%"

        if neurodegen_flag and "flybase_gene_id" in meta_df.columns:
            meta_df["GENE_HAS_NEURODEGEN_PAPER_REF"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_neurodegen_genes(fbgnid)
            )
            meta_df["GENE_NEURODEGEN_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_nd_paper_refs(fbgnid)
            )
            pct_nd = (meta_df[meta_df["GENE_HAS_NEURODEGEN_PAPER_REF"]]["flybase_gene_id"].nunique() / meta_df["flybase_gene_id"].nunique() * 100)
            print(f"Percentage of genes with neurodegeneration references: {pct_nd:.2f}%")
            summary_dict["% of genes with neurodegen paper ref"] = f"{pct_nd:.2f}%"

        if args.neurodegen_specific and "flybase_gene_id" in meta_df.columns:
            meta_df["GENE_HAS_NEURODEGEN_SPECIFIC_REF"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_neurodegen_specific_genes(fbgnid)
            )
            meta_df["GENE_NEURODEGEN_SPECIFIC_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_nd_specific_paper_refs(fbgnid)
            )
            pct_nd_specific = (
                    meta_df[meta_df["GENE_HAS_NEURODEGEN_SPECIFIC_REF"]]["flybase_gene_id"].nunique()
                    / meta_df["flybase_gene_id"].nunique() * 100
            )
            print(f"Percentage of genes with neurodegeneration-specific references: {pct_nd_specific:.2f}%")
            summary_dict["% of genes with neurodegen-specific paper ref"] = f"{pct_nd_specific:.2f}%"

        if args.TDP43 and "flybase_gene_id" in meta_df.columns:
            meta_df["GENE_HAS_TDP43_PAPER_REF"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_tdp43_genes(fbgnid)
            )
            meta_df["GENE_TDP43_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_tdp43_paper_refs(fbgnid)
            )
            pct_tdp43 = (
                    meta_df[meta_df["GENE_HAS_TDP43_PAPER_REF"]]["flybase_gene_id"].nunique()
                    / meta_df["flybase_gene_id"].nunique() * 100
            )
            print(f"Percentage of genes with TDP-43 references: {pct_tdp43:.2f}%")
            summary_dict["% of genes with TDP-43 paper ref"] = f"{pct_tdp43:.2f}%"

        summary_df_temp.append(summary_dict)


        # Reorder columns so the newly created columns end up at the right side
        new_cols = [
            "GENE_HAS_SLEEP_PAPER_REF", "GENE_SLEEP_PAPER_REFS",
            "GENE_HAS_NEURODEGEN_PAPER_REF", "GENE_NEURODEGEN_PAPER_REFS",
            "GENE_HAS_NEURODEGEN_SPECIFIC_REF", "GENE_NEURODEGEN_SPECIFIC_PAPER_REFS",
            "GENE_HAS_TDP43_PAPER_REF", "GENE_TDP43_PAPER_REFS"
        ]
        cols_to_order = [c for c in new_cols if c in meta_df.columns]
        other_cols = [c for c in meta_df.columns if c not in cols_to_order]
        meta_df = meta_df[other_cols + cols_to_order]

        # Check blacklists
        all_refs_cols = [col for col in cols_to_order if "PAPER_REFS" in col]
        for refs_col in all_refs_cols:
            all_refs = set()
            for val in meta_df[refs_col]:
                if isinstance(val, str) and val.strip() != '-':
                    for ref in val.split(','):
                        all_refs.add(ref.strip())
            intersection = all_refs.intersection(blacklisted_refs_set)
            if intersection:
                print(f"\n[ERROR] Found blacklisted references in column '{refs_col}':")
                print(intersection)
                print("Exiting script due to blacklisted references.")
                exit(1)

        def count_refs(ref_str):
            if not isinstance(ref_str, str) or ref_str.strip() == '-' or ref_str.strip() == '':
                return 0
            # Split by comma, ignoring extra spaces
            return len([x for x in ref_str.split(',') if x.strip() != ''])
        # Build DOIs & minirefs in newline format & fetch Crossref titles
        if "GENE_SLEEP_PAPER_REFS" in meta_df.columns:
            meta_df[["GENE_SLEEP_PAPER_DOI", "GENE_SLEEP_PAPER_MINIREF"]] = (
                meta_df["GENE_SLEEP_PAPER_REFS"].apply(
                    lambda ref_string: pd.Series(get_ref_info(ref_string, fbrf_paper_df))
                )
            )
            meta_df["GENE_SLEEP_PAPER_TITLE"] = meta_df["GENE_SLEEP_PAPER_DOI"].apply(fetch_titles_for_dois)

            meta_df['SLEEP_PAPER_REFS_COUNT'] = meta_df['GENE_SLEEP_PAPER_REFS'].apply(count_refs)

        if "GENE_NEURODEGEN_PAPER_REFS" in meta_df.columns:
            meta_df[["GENE_NEURODEGEN_PAPER_DOI", "GENE_NEURODEGEN_PAPER_MINIREF"]] = (
                meta_df["GENE_NEURODEGEN_PAPER_REFS"].apply(
                    lambda ref_string: pd.Series(get_ref_info(ref_string, fbrf_paper_df))
                )
            )
            meta_df["GENE_NEURODEGEN_PAPER_TITLE"] = meta_df["GENE_NEURODEGEN_PAPER_DOI"].apply(fetch_titles_for_dois)
            meta_df['ND_PAPER_REFS_COUNT'] = meta_df['GENE_NEURODEGEN_PAPER_REFS'].apply(count_refs)

        if "GENE_NEURODEGEN_SPECIFIC_PAPER_REFS" in meta_df.columns:
            meta_df[["GENE_NEURODEGEN_SPECIFIC_PAPER_DOI", "GENE_NEURODEGEN_SPECIFIC_PAPER_MINIREF"]] = (
                meta_df["GENE_NEURODEGEN_SPECIFIC_PAPER_REFS"].apply(
                    lambda ref_string: pd.Series(get_ref_info(ref_string, fbrf_paper_df))
                )
            )
            meta_df["GENE_NEURODEGEN_SPECIFIC_PAPER_TITLE"] = meta_df["GENE_NEURODEGEN_SPECIFIC_PAPER_DOI"].apply(fetch_titles_for_dois)
            meta_df['ND_SPECIFIC_PAPER_REFS_COUNT'] = meta_df['GENE_NEURODEGEN_SPECIFIC_PAPER_REFS'].apply(count_refs)

        if "GENE_HAS_TDP43_PAPER_REF" in meta_df.columns:
            meta_df[["GENE_TDP43_PAPER_DOI", "GENE_TDP43_PAPER_MINIREF"]] = (
                meta_df["GENE_TDP43_PAPER_REFS"].apply(
                    lambda ref_string: pd.Series(get_ref_info(ref_string, fbrf_paper_df))
                )
            )
            meta_df["GENE_TDP43_PAPER_TITLE"] = meta_df["GENE_TDP43_PAPER_DOI"].apply(fetch_titles_for_dois)
            meta_df['TDP43_PAPER_REFS_COUNT'] = meta_df['GENE_TDP43_PAPER_REFS'].apply(count_refs)

        # Sort if needed
        if sortby and sortby in meta_df.columns:
            meta_df[sortby] = pd.to_numeric(meta_df[sortby], errors='coerce').fillna(0)
            meta_df.sort_values(sortby, ascending=False, inplace=True)

        meta_df.fillna("-", inplace=True)
        meta_df = meta_df.loc[:, ~meta_df.columns.str.contains('Unnamed')]
        # regex = r'(DGE_|DEG_|Rain_)'
        # meta_df = meta_df.loc[:, ~meta_df.columns.str.contains(regex, case=False)]
        meta_df.to_csv(target_file, index=False)
    # Finally, create summary CSV
    summary_df = pd.DataFrame(summary_df_temp)
    if not summary_df.empty:
        summary_df.sort_values("Total Genes", ascending=False, inplace=True)
        summary_df.to_csv(os.path.join(target_dir, "gene_paper_refs_summary_table.csv"), index=False)
    cache_items = list(doi_title_cache.items())  # [(doi, title), (doi, title), ...]
    df_to_save = pd.DataFrame(cache_items, columns=["doi", "title"])
    df_to_save.to_csv(CACHE_FILE, index=False)
    print(f"[INFO] Saved updated cache of {len(doi_title_cache)} DOIs to {CACHE_FILE}.")
