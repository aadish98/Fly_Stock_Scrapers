import pandas as pd
import pandas as pd
import sys
import os
import glob
from suds.client import Client
import csv
import numpy as np

import warnings
warnings.filterwarnings("ignore")

# Read all FlyBase published precompiled tsv files
folder_path = 'Data/FlyBase/'
tsv_files = glob.glob(os.path.join(folder_path, '**/*.tsv'), recursive=True)
flybase_dict = {}
for file in tsv_files:
    file_name = os.path.splitext(os.path.basename(file))[0]
    df = pd.read_csv(file, sep='\t', on_bad_lines='skip', low_memory= False)
    flybase_dict[file_name] = df
fbrf_paper_set = set(flybase_dict["fbrf_pmid_pmcid_doi_fb_2024_04"].loc[
    flybase_dict["fbrf_pmid_pmcid_doi_fb_2024_04"]['pub_type'] == 'paper', 'FBrf'])

# Read the file (assuming it's in a tab-separated format)

directory_path = sys.argv[1]
if not os.path.exists(directory_path):
    print("ERROR: Folder path entered does not exist!")

csv_files = glob.glob(os.path.join(directory_path, '**.csv'), recursive=True)
summary_list = []
row_dict = {}
# Define custom logic for processing <newline> as actual newlines within fields
counter = 0
for csv_file in csv_files:
    row_dict = {}
    if("summary" in csv_file):
        continue
    print(f'Processing file {csv_file}')
    df = pd.read_csv(csv_file)[["RefIDs"]].drop_duplicates(subset= "RefIDs")
    print("Unique input: ")
    print(df.shape[0])
    df = df.merge(flybase_dict["entity_publication_fb_2024_04"][["FlyBase_publication_id", "entity_id", "entity_name"]], left_on="RefIDs", right_on = "FlyBase_publication_id")
    #print(df.columns)
    df.drop(columns= ['FlyBase_publication_id'], inplace = True)
    print("All entity id's: ")
    print(df.shape[0])
    df = df.drop_duplicates(subset= "entity_name")
    print("Unique id's: ")
    print(df.shape[0])
    #df = df.merge(flybase_dict["Alleles_2_Gene_fbal_to_fbgn_fb_2024_04"], left_on = "entity_name", right_on = "AlleleSymbol")
    #df.drop_duplicates(subset="AlleleID", inplace = True)
    print(df.shape[0])
    df['category'] = df['entity_id'].str.extract(r'^(FB\w{2})')

# Find unique categories
    unique_categories = df['category'].unique()

    # Display the unique categories
    print(unique_categories)
    df.to_csv(csv_file, index=True)