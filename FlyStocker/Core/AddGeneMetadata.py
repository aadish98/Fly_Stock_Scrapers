# Example Usage:
# $ python3 AddStockMetadata.py <InputDirectory> <OutputDirectory>
# Function: Adds experiment results/ better refs in processed csv (derived from flybase batchdownload tsv)

import pandas as pd
import os
import glob
import numpy as np
import sys
import argparse

blacklisted_refs = pd.read_csv("../../Data/FlyBase/custom/blacklisted_refs.csv")
blacklisted_refs_set = set(blacklisted_refs["FBrf"])
# ----
fbrf_sleep_genes_set = pd.read_csv(
    "../../Data/FlyBase/QueryResults/Cross_Ref_Tables/Sleep/sleepPaperQuery_crossRef_table.csv"
)
fbrf_sleep_genes_set["blacklisted_ref"] = fbrf_sleep_genes_set["ref_id"].apply(
    lambda x: True if x in set(blacklisted_refs["FBrf"]) else False
)
fbrf_sleep_genes_set = fbrf_sleep_genes_set[fbrf_sleep_genes_set["blacklisted_ref"] == False]
# ---
fbrf_neurodegen_genes_set = pd.read_csv(
    "../../Data/FlyBase/QueryResults/Cross_Ref_Tables/Neurodegen/Neurodegen-Broad/neurodegenPaperQuery_crossRef_table.csv"
)
fbrf_neurodegen_genes_set["blacklisted_ref"] = fbrf_neurodegen_genes_set["ref_id"].apply(
    lambda x: True if x in set(blacklisted_refs["FBrf"]) else False
)
fbrf_neurodegen_genes_set = fbrf_neurodegen_genes_set[fbrf_neurodegen_genes_set["blacklisted_ref"] == False]
# ---
fbrf_neurodegen_specific_genes_set = pd.read_csv(
    "../../Data/FlyBase/QueryResults/Cross_Ref_Tables/Neurodegen/Neurodegen-Specific/neurodegenPaperQuery_crossRef_table.csv"
)
fbrf_neurodegen_specific_genes_set["blacklisted_ref"] = fbrf_neurodegen_specific_genes_set["ref_id"].apply(
    lambda x: True if x in set(blacklisted_refs["FBrf"]) else False
)
fbrf_neurodegen_specific_genes_set = fbrf_neurodegen_specific_genes_set[fbrf_neurodegen_specific_genes_set["blacklisted_ref"] == False]
# ---
fbrf_tdp43_genes_set = pd.read_csv(
    "../../Data/FlyBase/QueryResults/Cross_Ref_Tables/Neurodegen/TDP-43/TDP-43PaperQuery_crossRef_table.csv"
)
fbrf_tdp43_genes_set["blacklisted_ref"] = fbrf_tdp43_genes_set["ref_id"].apply(
    lambda x: True if x in set(blacklisted_refs["FBrf"]) else False
)
fbrf_tdp43_genes_set = fbrf_tdp43_genes_set[fbrf_tdp43_genes_set["blacklisted_ref"] == False]

def remove_blacklisted_refs(ref_string, blacklist):
    """
    Takes a comma-separated string of references, removes any that are in `blacklist`,
    and returns a cleaned, comma-separated string (or NaN if empty).
    """
    if pd.isna(ref_string):
        return np.nan  # nothing to clean
    refs = [r.strip() for r in ref_string.split(",")]
    filtered = [r for r in refs if r not in blacklist and r != ""]
    if not filtered:
        return np.nan  # if nothing left, store as NaN
    return ", ".join(filtered)

# ---
def get_sleep_genes(fbgnid):
    sleep_genes_set = set(fbrf_sleep_genes_set["flybase_gene_id"].astype(str))  # Ensure set has str values
    if not isinstance(fbgnid, str) or pd.isnull(fbgnid) or fbgnid == '-':
        return False
    if fbgnid in sleep_genes_set:
        return True
    else:
        return False

def get_sleep_paper_refs(fbgnid):
    sleep_refs = fbrf_sleep_genes_set[fbrf_sleep_genes_set["flybase_gene_id"]== fbgnid]["ref_id"].unique()
    sleep_refs_str = ', '.join(sleep_refs.astype(str).tolist())
    return sleep_refs_str if sleep_refs_str else '-'

# ---
def get_neurodegen_genes(fbgnid):
    nd_genes_set = set(fbrf_neurodegen_genes_set["flybase_gene_id"].astype(str))
    if not isinstance(fbgnid, str) or pd.isnull(fbgnid) or fbgnid == '-':
        return False
    if fbgnid in nd_genes_set:
        return True
    else:
        return False

def get_nd_paper_refs(fbgnid):
    nd_refs = fbrf_neurodegen_genes_set[fbrf_neurodegen_genes_set["flybase_gene_id"]== fbgnid]["ref_id"].unique()
    nd_refs_str = ', '.join(nd_refs.astype(str).tolist())
    return nd_refs_str if nd_refs_str else '-'
# ---
def get_neurodegen_specific_genes(fbgnid):
    nd_genes_set = set(fbrf_neurodegen_specific_genes_set["flybase_gene_id"].astype(str))
    if not isinstance(fbgnid, str) or pd.isnull(fbgnid) or fbgnid == '-':
        return False
    if fbgnid in nd_genes_set:
        return True
    else:
        return False

def get_nd_specific_paper_refs(fbgnid):
    nd_refs = fbrf_neurodegen_specific_genes_set[fbrf_neurodegen_specific_genes_set["flybase_gene_id"]== fbgnid]["ref_id"].unique()
    nd_refs_str = ', '.join(nd_refs.astype(str).tolist())
    return nd_refs_str if nd_refs_str else '-'
# ---
def get_tdp43_genes(fbgnid):
    nd_genes_set = set(fbrf_tdp43_genes_set["flybase_gene_id"].astype(str))
    if not isinstance(fbgnid, str) or pd.isnull(fbgnid) or fbgnid == '-':
        return False
    if fbgnid in nd_genes_set:
        return True
    else:
        return False

def get_tdp43_paper_refs(fbgnid):
    nd_refs = fbrf_tdp43_genes_set[fbrf_tdp43_genes_set["flybase_gene_id"]== fbgnid]["ref_id"].unique()
    nd_refs_str = ', '.join(nd_refs.astype(str).tolist())
    return nd_refs_str if nd_refs_str else '-'

if __name__ == "__main__":
    # Get directory paths and optional gene column name from command-line arguments
    # Reference experiment results files (consistent sleep/ wake association)
    parser = argparse.ArgumentParser()
    parser.add_argument("directory_path")
    parser.add_argument("target_dir")
    parser.add_argument("mergeOn")
    parser.add_argument("-neurodegen", action="store_true", help="Set to true if neurodegeneration is required")
    parser.add_argument("-sleep", action="store_true", help="Set to true if neurodegeneration is required")
    parser.add_argument("-sortby", choices=["frequency"], help="Specify sorting criteria")
    args = parser.parse_args()

    directory_path = args.directory_path
    target_dir = args.target_dir
    neurodegen_flag = args.neurodegen
    sleep_flag = args.sleep
    sortby = args.sortby
    mergeOn = args.mergeOn
    csv_files = glob.glob(os.path.join(directory_path, "**.csv"), recursive=True)
    print(csv_files)
    df_list = []

    print("Reference Files: \n")
    # Loop through the list of CSV files and read each into a DataFrame
    for file in csv_files:
        print(repr(file))
        print('\n')
        df = pd.read_csv(file, keep_default_na=False)  # Read each CSV file
        df_list.append(df)  # Append the DataFrame to the list
    # Concatenate all the DataFrames in the list into a single DataFrame
    combined_df = pd.concat(df_list, ignore_index=True)
    duplicate_counts = combined_df["flybase_gene_id"].value_counts()
    if "flybase_gene_id" not in combined_df.columns:
        raise RuntimeError("METADATA DF DOESNT CONTAIN FBgnID!!!")
    # if "AlleleID" in combined_df.columns:
    #     combined_df.drop_duplicates(subset="AlleleID", inplace=True)
    target_files = glob.glob(os.path.join(target_dir, "**.csv"), recursive=True)
    print(target_files)

    #allele_id_to_refs = refs_entity_mapping.groupby("entity_id")["FlyBase_publication_id"].apply(set).to_dict()
    for target_file in target_files:
        target_df = pd.read_csv(target_file, keep_default_na=False)

        # 2) Preprocess the stock gene metadata- in this case experiment results- such that cols arent duplicated
        # (i) Remove all duplicated cols except AlleleID since this is our matching key
        common_columns = list(
            set(target_df.columns).intersection(set(combined_df.columns))
        )
        if "AlleleID" in common_columns:
            common_columns.remove("AlleleID")
        if "flybase_gene_id" in common_columns:
            common_columns.remove("flybase_gene_id")
        # (ii) Add any trailing cols with Unnamed in them
        common_columns += [
            col for col in combined_df.columns if "unnamed" in col.lower()
        ]
        combined_df.drop(columns=common_columns, inplace=True)
        if "flybase_gene_id_y" in target_df.columns:
            target_df["flybase_gene_id"] = target_df["flybase_gene_id_x"]
            target_df.drop(columns=["flybase_gene_id_y", "flybase_gene_id_x"], inplace = True)

        if "AlleleID_x" in target_df.columns:
            target_df.drop(columns=["AlleleID_x", "AlleleID_y"], inplace = True)
        print(f"Original file length: {target_df.shape[0]}")

        ## CHANGE FOR ADDING METADATA AS PER ALLELEID OR FLYBASE ID
        meta_df = pd.merge(target_df, combined_df, how="left", on = mergeOn)
        if "flybase_gene_id_y" in meta_df.columns:
            meta_df["flybase_gene_id"] = meta_df["flybase_gene_id_x"]
            meta_df.drop(columns=["flybase_gene_id_y", "flybase_gene_id_x"], inplace = True)

        if "AlleleID_x" in meta_df.columns:
            meta_df.drop(columns=["AlleleID_x", "AlleleID_y"], inplace = True)
        # 3i) Check for gene in relevant published gene sets
        if sleep_flag:
            meta_df["GENE_HAS_SLEEP_PAPER_REF"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_sleep_genes(fbgnid)
            )
            meta_df["GENE_SLEEP_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_sleep_paper_refs(fbgnid)
            )
            meta_df.loc[meta_df["GENE_SLEEP_PAPER_REFS"] == "-", "GENE_HAS_SLEEP_PAPER_REF"] = False
        if neurodegen_flag:
            meta_df["GENE_HAS_NEURODEGEN_PAPER_REF"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_neurodegen_genes(fbgnid)
            )
            meta_df["GENE_NEURODEGEN_PAPER_REFS"] = meta_df["flybase_gene_id"].apply(
                lambda fbgnid: get_nd_paper_refs(fbgnid)
            )

            meta_df.loc[meta_df["GENE_NEURODEGEN_PAPER_REFS"] == "-", "GENE_HAS_NEURODEGEN_PAPER_REF"] = False
        # 3ii) Remove all blacklisted refs from the 'GENE refs' column
        cols = [
            ("GENE_HAS_SLEEP_PAPER_REF", "GENE_SLEEP_PAPER_REFS"),
            ("GENE_HAS_NEURODEGEN_PAPER_REF", "GENE_NEURODEGEN_PAPER_REFS")
        ]

        for has_col, refs_col in cols:
            # Only proceed if both columns actually exist in df
            if has_col in meta_df.columns and refs_col in meta_df.columns:
                # 1. Clean the references by removing blacklisted items
                meta_df[refs_col] = meta_df[refs_col].apply(remove_blacklisted_refs, blacklist=blacklisted_refs_set)
                # 2. Recompute the boolean: True if REFS column is not empty/NaN
                meta_df[has_col] = meta_df[refs_col].notna()
        cols_to_order = [col for col in ["STOCK_PAPER_REFS", "GENE_HAS_SLEEP_PAPER_REF","GENE_SLEEP_PAPER_REFS", "GENE_HAS_NEURODEGEN_PAPER_REF", "GENE_NEURODEGEN_PAPER_REFS"] if col in meta_df.columns]

        if sleep_flag:
            meta_df.loc[meta_df["GENE_SLEEP_PAPER_REFS"] == "-", "GENE_HAS_SLEEP_PAPER_REF"] = False
        if neurodegen_flag:
            meta_df.loc[meta_df["GENE_NEURODEGEN_PAPER_REFS"] == "-", "GENE_HAS_NEURODEGEN_PAPER_REF"] = False
        # Reorder the columns
        other_cols = [col for col in meta_df.columns if col not in cols_to_order]
        meta_df = meta_df[other_cols + cols_to_order]
        # Remove duplicate rows based on 'StockID' and 'collection_short_name' columns
        meta_df = meta_df.drop_duplicates(subset=['FBst'])

        # 4) Sort by frequency if specified
        if sortby == "frequency":
            sorted_df = meta_df.sort_values(by="frequency", ascending=False)
            sorted_df.to_csv(target_file, index=False)
        else:
            meta_df.to_csv(target_file, index=True)

        print(f"Final length: {meta_df.shape[0]}")
        print(
            f'Number of stocks with no metadata: {len(meta_df[pd.isnull(meta_df["flybase_gene_id"])])}'
        )
