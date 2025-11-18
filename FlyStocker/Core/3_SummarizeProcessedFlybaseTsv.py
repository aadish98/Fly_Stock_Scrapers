# Usage: python3 3_SummarizeProcessedFlybaseTsv.py /path/to/process_csv_files /path/to/config.json
import pandas as pd
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import warnings
import argparse
import re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from HelperScripts.plotting_functions import *
from HelperScripts.processing_helper_functions import *
from collections import defaultdict

warnings.filterwarnings("ignore")
plt.rcParams["figure.autolayout"] = True

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Process stock data CSV files with configurable filters."
    )
    parser.add_argument(
        "directory_path", type=str, help="Path to the directory containing CSV files"
    )
    parser.add_argument(
        "config_path", type=str, help="Path to the configuration JSON file"
    )

    parser.add_argument(
        "--relevantSearchTerm",
        type=str,
        help="sleep or neurodegen",
    )
    parser.add_argument(
        "--maxStocksPerGene",
        type=int,
        default=5,
        help="Maximum number of stocks to include per gene (default: 5)",
    )
    parser.add_argument(
        "--maxStocksPerAllele",
        type=int,
        default=3,
        help="Maximum number of stocks to include per allele (default: 3)",
    )
    args = parser.parse_args()

    directory_path = args.directory_path
    config = load_split_config(args.config_path)
    relevantSearchTerm = args.relevantSearchTerm
    maxStocksPerGene = args.maxStocksPerGene
    maxStocksPerAllele = args.maxStocksPerAllele

    if not os.path.exists(directory_path):
        print("ERROR: Folder path entered does not exist!")
    flybase_genes_path = "../../Data/FlyBase/Genes/fb_synonym_fb_2025_03.tsv"
    flybase_genes = pd.read_csv(flybase_genes_path, sep="\t")

    bdsc_stock_table_fp = "../../Data/BDSC/bloomington.csv"
    bdsc_stocks = pd.read_csv(
        bdsc_stock_table_fp, usecols=["Stk #", "All Affected Chromosomes"]
    )

    # Prepare a set of all possible gene symbols and their synonyms
    flybase_genes = flybase_genes[
        flybase_genes["primary_FBid"].str.contains("FBgn", na=False)
    ]
    flybase_genes["symbol_synonym_list"] = flybase_genes["symbol_synonym(s)"].str.split("|")

    # Initialize an empty set for all gene symbols and synonyms
    all_gene_symbols = set()

    # Populate the set with current symbols and synonyms
    for _, row in flybase_genes.iterrows():
        # Add the current symbol
        all_gene_symbols.add(row["current_symbol"])
        # Add all synonyms
        if isinstance(row["symbol_synonym_list"], list):
            all_gene_symbols.update(row["symbol_synonym_list"])

    stockgenes = pd.read_csv("../../Data/BDSC/stockgenes.csv")
    stockcomps_map_comments = pd.read_csv("../../Data/BDSC/stockcomps_map_comments.csv")
    stockcomps_map_comments = stockcomps_map_comments.rename(columns= {'Genotype': "BDSC_Genotype"})
    # Set up command-line argument parsing

    filters_config = config["filters"]
    combinations = config["combinations"]

    balancers_file_path = "../../Data/Metadata/balancers.csv"
    balancers_df = pd.read_csv(balancers_file_path)
    balancer_names = balancers_df["balancer_name"].tolist()
    insertion_mapping_df = read_insertion_mapping_df()
    gene_mapping_df = read_gene_mapping_df()
    csv_files = glob.glob(os.path.join(args.directory_path, "*.csv"), recursive=True)
    # Define custom logic for processing <newline> as actual newlines within fields


    overall_summary_dict = defaultdict(dict)
    counter = 0
    for csv_file in csv_files:

        if "summary" in csv_file or "enriched_stocker_file" in csv_file:
            continue
        # Create the necessary folder for the breakdown
        if not os.path.exists(os.path.join(directory_path, os.path.basename(csv_file))):
            os.makedirs(os.path.join(directory_path, os.path.basename(csv_file)))
        stock_sheet_folder = os.path.join(
            directory_path, os.path.basename(csv_file).replace(".csv", "")
        )
        if not os.path.exists(os.path.join(stock_sheet_folder, "Stock_Breakdown")):
            os.makedirs(
                os.path.join(os.path.join(stock_sheet_folder, "Stock_Breakdown")),
                exist_ok=True,
            )
        if not os.path.exists(
            os.path.join(stock_sheet_folder, "All_Stocks_Summary_Plots")
        ):
            os.makedirs(
                os.path.join(
                    os.path.join(stock_sheet_folder, "All_Stocks_Summary_Plots")
                ),
                exist_ok=True,
            )
        print(f"Processing file {csv_file}")

        df = pd.read_csv(csv_file)
        if 'GO_Terms' in df.columns:
            sleep_GO_terms = set(pd.read_csv(
                "/Users/AadishShah/Desktop/Projects.nosync/Allada Lab Downloads/Fly_Stock_Scrapers/Data/Metadata/sleep_GO_Terms.csv")[
                                     "Term"])
            wake_GO_terms = set(pd.read_csv(
                "/Users/AadishShah/Desktop/Projects.nosync/Allada Lab Downloads/Fly_Stock_Scrapers/Data/Metadata/wake_GO_Terms.csv")[
                                    "Term"])
        df = df[(df["StockID"] != "-") & (df["StockID"] != "")]
        num_unique_genes = len(np.unique(df["flybase_gene_id"].astype(str)))
        unique_flybase_gene_ids = np.unique(df["flybase_gene_id"])
        fbgnid_counts_dict = {}
        for fbgnid in unique_flybase_gene_ids:
            fbgnid_counts_dict[fbgnid] = 0
        unique_allele_ids = np.unique(df["flybase_allele_id"])
        alleleid_counts_dict = {}
        for allele_id in unique_allele_ids:
            alleleid_counts_dict[allele_id] = 0
        stockid_set = set()
        if "STOCK_PAPER_REFS" in df.columns:
            df = df.rename({"STOCK_PAPER_REFS": "ALLELE_PAPER_REFS"})
        df = drop_suffix_columns(df)
        df = df.drop(
            columns=list(
                set(
                    insertion_mapping_df.columns.tolist()
                    + gene_mapping_df.columns.tolist()
                )
            ),
            errors="ignore",
        )
        df["AlleleID"] = df["flybase_allele_id"]
        df["Balancers"] = df["FB_genotype"].apply(
            lambda genotype: ", ".join(
                [balancer for balancer in balancer_names if balancer in genotype]
            )
            or "-"
        )

        df["All Affected Chromosomes"] = df.apply(
            lambda row: (
                affected_chromosomes(row["StockID"], bdsc_stocks)
                if "Bloomington" in row["collection_short_name"]
                else ""
            ),
            axis=1,
        )

        # Add insertion specific information (chromosome arm/ location)
        df["insertion_id"] = df["ASSOCIATED_INSERTION"].apply(
            lambda x: (
                re.search(r"^<fbsym>(FBti\d+)", x).group(1)
                if re.search(r"^<fbsym>(FBti\d+)", x)
                else x
            )
        )
        df["insertion_id"] = df["insertion_id"].apply(
            lambda x: x.strip() if isinstance(x, str) else x
        )
        df = df.drop(
            columns=list(
                set(
                    ["genomic_location", "chromosome_arm", "position_range"]
                    + insertion_mapping_df.columns.tolist()
                )
            ),
            errors="ignore",
        )
        df = insertion_join(df, insertion_mapping_df)
        df = df.rename(
            columns={"genomic_location": "ASSOCIATED_INSERTION_genomic_location"}
        )
        df[["insertion_chromosome_arm", "insertion_position_range"]] = df[
            "ASSOCIATED_INSERTION_genomic_location"
        ].str.split(":", expand=True)
        # Insertion specific genomic location data
        df = df.drop(
            columns=[
                "range",
                "orientation",
                "estimated_cytogenetic_location",
                "observed_cytogenetic_location",
                "join_key",
            ],
            errors="ignore",
        )  # Comment to keep these cols
        # Add gene specific location data

        df = df.merge(
            gene_mapping_df, left_on="flybase_gene_id", right_on="primary_FBid"
        )
        df = df.drop(columns=["primary_FBid"])
        df["recombination_loc"] = df["recombination_loc"].astype(str)
        df = df.rename(
            columns={
                "sequence_loc": "gene_sequence_loc",
                "recombination_loc": "gene_recombination_loc",
                "cytogenetic_loc": "gene_cytogenetic_loc",
            }
        )
        df["gene_recombination_loc"] = df["gene_recombination_loc"].astype(str)
        df["gene_recombination_loc"] = df["gene_recombination_loc"].apply(
            lambda x: f"'{x}'" if isinstance(x, str) else x
        )


        df["gene_cytogenetic_loc"] = df["gene_cytogenetic_loc"].astype(str)
        df[["gene_chromosome_arm", "gene_location"]] = df[
            "gene_sequence_loc"
        ].str.split(":", expand=True)
        df = df.drop(columns=gene_mapping_df.columns.tolist(), errors="ignore")
        if "GENE_SLEEP_PAPER_REFS" in df.columns:
            df["GENE_SLEEP_PAPER_REFS"] = df["GENE_SLEEP_PAPER_REFS"].astype(str)
        if relevantSearchTerm == 'sleep':
            df["ALLELE_SLEEP_PAPERS"] = df.apply(
                lambda row: find_overlaps(
                    row["ALLELE_PAPER_REFS"], row["GENE_SLEEP_PAPER_REFS"]
                ),
                axis=1,
            )

        if relevantSearchTerm == 'neurodegen':
            df["ALLELE_NEURODEGEN_PAPERS"] = df.apply(
                lambda row: find_overlaps(
                    row["ALLELE_PAPER_REFS"], row["GENE_NEURODEGEN_PAPER_REFS"]
                ),
                axis=1,
            )

        relevance_columns = {
            "sleep": "ALLELE_SLEEP_PAPERS",
            "neurodegen": "ALLELE_NEURODEGEN_PAPERS",
        }
        if relevantSearchTerm in relevance_columns.keys():
            df["ALLELE_PAPER_RELEVANCE_SCORE"] = df.apply(
                lambda row: calculate_relevance_score(
                    row[relevance_columns[relevantSearchTerm]], row["ALLELE_PAPER_REFS"]
                ),
                axis=1,
            )

        df["NUM_ALLELE_PAPER_REFS"] = df["ALLELE_PAPER_REFS"].apply(
            lambda x: 0 if x == "" else len(x.split(", "))
        )
        print("OG")
        print(df.shape)
        print(df["FBst"].nunique())
        df["StockID"] = df["StockID"].astype(str)
        stockcomps_map_comments['Stk #'] = stockcomps_map_comments['Stk #'].astype(str)
        df_bloom = df[df["collection_short_name"] == "Bloomington"]
        df_non_bloom = df[df["collection_short_name"] != "Bloomington"]
        print(df_bloom["FBst"].nunique() + df_non_bloom["FBst"].nunique())
        # Perform the merge only on the Bloomington subset


        df_bloom_merged = df_bloom.merge(
            stockcomps_map_comments[["Stk #", "comment1", "BDSC_Genotype"]],
            left_on=["StockID"],
            right_on=["Stk #"],
            how="left",
        )
        duplicated_stockids = df_bloom_merged[df_bloom_merged["StockID"].duplicated(keep=False)]

        # Recombine the Bloomington-merged subset with the non-Bloomington subset
        df = pd.concat([df_bloom_merged, df_non_bloom], ignore_index=True)
        print("After first merge")
        print(df.shape)
        print(df["FBst"].nunique())
        df = df.sort_values(by="NUM_ALLELE_PAPER_REFS", ascending=False)

        # Drop suffix cols i.e. _x or _y
        df = drop_suffix_columns(df)
        print(df.shape)
        df["multiple_insertions"] = df["Genotype"].str.count("attP") > 1
        if 'GO_Terms' in df.columns:
            if 'sleep' in csv_file:
                df[['GO_plus', 'Relevant_GO_Terms']] = df['GO_Terms'].apply(lambda x: find_relevant_go_terms(x, sleep_GO_terms))
            if 'wake' in csv_file:
                df[['GO_plus', 'Relevant_GO_Terms']] = df['GO_Terms'].apply(lambda x: find_relevant_go_terms(x, wake_GO_terms))
        new_file_name = csv_file.replace(".csv", "_enriched_stocker_file.csv")
        stock_subfolders = {
            "all_stocks": os.path.join(stock_sheet_folder, "ALL_Stocks"),
            "special_cases": os.path.join(stock_sheet_folder, "Special_Cases"),
        }
        os.makedirs(stock_subfolders["all_stocks"], exist_ok=True)

        new_file_path = os.path.join(
            stock_subfolders["all_stocks"], os.path.basename(new_file_name)
        )
        df = df.loc[:, ~df.columns.str.contains(r'^Unnamed')]
        df.to_csv(new_file_path, index=False)  # Save all stocks found by algo

        df = df.fillna('-')
        df = df.drop_duplicates()
        filters = apply_filters(df, filters_config)
        # Initialize an empty list to store the summary data
        summary_data = []
        # Loop over each combination, filter the dataset, and generate the summary
        base_folder = os.path.join(stock_sheet_folder, "Stock_Breakdown")
        stock_paper_relevance_map = {
            0: "stock has no paper ref",
            1: "stock has paper ref",
            2: f'stock has {relevantSearchTerm} paper ref',
        }

        for combination in combinations:
            filtered_stocks = df
            for layer in combination:
                filter_config = filters_config[layer]
                column = filter_config["column"]
                filter_type = filter_config["type"]
                value = filter_config["value"]
                if filter_type == "contains":
                    filtered_stocks = filtered_stocks[filtered_stocks[column].str.contains(value, na=False)]
                elif filter_type == "doesnt_contain":
                    filtered_stocks = filtered_stocks[~filtered_stocks[column].str.contains(value, na=False)]
                elif filter_type == "equals":
                    filtered_stocks = filtered_stocks[filtered_stocks[column] == value]
                elif filter_type == "doesnt_equal":
                    filtered_stocks = filtered_stocks[filtered_stocks[column] != value]
                elif filter_type == "insertion_site":
                    filtered_stocks = filtered_stocks[(filtered_stocks[column].str.contains(value, case = True, na = False)) & (~filtered_stocks["multiple_insertions"])]
            if "join_key" in filtered_stocks.columns:
                filtered_stocks = filtered_stocks.drop(
                    columns=["join_key"], errors="ignore"
                )
            filtered_stocks = filtered_stocks.drop_duplicates(subset = ["FBst"])
            if "ALLELE_PAPER_RELEVANCE_SCORE" in filtered_stocks.columns:
                filtered_stocks = filtered_stocks.sort_values(
                        by='ALLELE_PAPER_RELEVANCE_SCORE',
                        ascending=False
                )
            if relevantSearchTerm == "sleep":
                filtered_stocks = filtered_stocks.sort_values(
                    by='ALLELE_SLEEP_PAPERS',
                    key=lambda col: col.apply(lambda x: 0 if x.strip() == '-' else len(x.split(','))),
                    ascending=False
                )
            if relevantSearchTerm == "neurodegen":
                filtered_stocks = filtered_stocks.sort_values(
                    by='ALLELE_NEURODEGEN_PAPERS',
                    key=lambda col: col.apply(lambda x: 0 if x.strip() == '-' else len(x.split(','))),
                    ascending=False
                )

            # Create subfolder structure based on combination elements
            subfolder_path = os.path.join(
                base_folder, *[layer.replace(" ", "_") for layer in combination]
            )

            # Ensure the directory exists
            os.makedirs(subfolder_path, exist_ok=True)

            # Save filtered stocks to CSV
            output_filename = f"{'_'.join(combination)}_stocks.csv".replace(" ", "")
            output_path = os.path.join(subfolder_path, output_filename)
            row_data = {
                "Category": " >> ".join(
                    [layer.replace(" ", "_") for layer in combination]
                )
            }
            # Write the filtered stocks end point
            limited_df = pd.DataFrame(columns=filtered_stocks.columns)
            for i in range(filtered_stocks.shape[0]):
                current_stock_id = filtered_stocks["FBst"].tolist()[i]
                ind_row = pd.DataFrame(filtered_stocks[filtered_stocks["FBst"] == current_stock_id])
                if ind_row.shape[0] > 1:
                    print(ind_row)
                current_fbgnid = ind_row["flybase_gene_id"].tolist()[0]
                current_alleleid = ind_row["flybase_allele_id"].tolist()[0]
                if (fbgnid_counts_dict[current_fbgnid] < maxStocksPerGene) and (alleleid_counts_dict[current_alleleid] < maxStocksPerAllele) and (current_stock_id not in stockid_set):
                    limited_df = pd.concat([limited_df, ind_row], ignore_index=True)
                    fbgnid_counts_dict[current_fbgnid] = fbgnid_counts_dict[current_fbgnid] + 1
                    alleleid_counts_dict[current_alleleid] = alleleid_counts_dict[current_alleleid] + 1
                    stockid_set.add(current_stock_id)

            if "frequency" in limited_df.columns:
                limited_df["frequency"] = limited_df["frequency"].replace("-", 0)
                limited_df["frequency"] = limited_df["frequency"].astype(float).astype(int)
                limited_df = limited_df.sort_values("frequency", ascending=False)
            else:
                print(limited_df)
                print(limited_df.shape[0])
            limited_df = limited_df.loc[:, ~limited_df.columns.str.contains(r'^Unnamed')]
            limited_df.to_csv(output_path)
            score_filtered = limited_df
            num_stocks = score_filtered["FBst"].nunique()
            num_unique_genes = score_filtered["flybase_gene_id"].nunique()
            num_unique_alleles = score_filtered["AlleleID"].nunique()
            row_data[
                f"Num Stocks | Num Alleles | Num Genes"
            ] = f"{num_stocks}  |  {num_unique_alleles} |  {num_unique_genes}"
            summary_data.append(row_data) # Adds combination summary
            if len(combination) >= 2:
                category_name = " >> ".join([combination[-2], combination[-1]])
            else:
                category_name = " >> ".join(combination)

            base_name = os.path.basename(csv_file)
            overall_summary_dict[base_name][f"Num_Stocks_{category_name}"] = num_stocks
            overall_summary_dict[base_name][f"Num_Alleles_{category_name}"] = num_unique_alleles
            overall_summary_dict[base_name][f"Num_Genes_{category_name}"] = num_unique_genes
        # Convert the summary data to a DataFrame
        summary_df = pd.DataFrame(summary_data)
        # Write summary table with three sheets based on STOCK_PAPER_RELEVANCE_SCORE
        summary_output_path = os.path.join(
            base_folder, "stock_files_summary_table.xlsx"
        )

        with pd.ExcelWriter(summary_output_path, engine="xlsxwriter") as writer:
            summary_df.to_excel(writer, index=False)
        # Add the histogram plots to the relevant dirs
        stock_sheet_folder = os.path.join(
            stock_sheet_folder, "All_Stocks_Summary_Plots"
        )
        plot_stock_distribution(df, csv_file, stock_sheet_folder)
        plot_stock_distribution(
            df, csv_file, stock_sheet_folder, fieldOfInterest="AlleleSymbol"
        )
    overall_summary_df = pd.DataFrame.from_dict(overall_summary_dict, orient="index").reset_index()
    overall_summary_df.rename(columns={"index": "File_Name"}, inplace=True)
    overall_summary_output_path = os.path.join(directory_path, "Overall_Stock_Summary_Table.xlsx")
    overall_summary_df.to_excel(overall_summary_output_path, index=False)

