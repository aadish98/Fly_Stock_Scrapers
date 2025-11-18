import pandas as pd
import sys
import os
import glob
import csv
import re
import warnings
import subprocess
# 1 cmdline arg: directory_path
warnings.filterwarnings("ignore")

csv.field_size_limit(sys.maxsize)

folder_path = '../../Data/FlyBase/'
tsv_files = glob.glob(os.path.join(folder_path, '**/*.tsv'), recursive=True)
flybase_dict = {}
for file in tsv_files:
    file_name = os.path.splitext(os.path.basename(file))[0]
    print(file_name)
    df = pd.read_csv(file, delimiter='\t', low_memory= False)
    flybase_dict[file_name] = df
fbrf_paper_set = set(flybase_dict["fbrf_pmid_pmcid_doi_fb_2024_06"].loc[
    flybase_dict["fbrf_pmid_pmcid_doi_fb_2024_06"]['pub_type'] == 'paper', 'FBrf'])

blacklisted_refs = pd.read_csv("../../Data/FlyBase/custom/blacklisted_refs.csv")
blacklisted_refs_set = set(blacklisted_refs["FBrf"])

def clean_references(refs):
    # Split by newlines and strip whitespace
    ref_list = []
    if '\n' in refs:
        ref_list = [ref.strip() for ref in refs.split("\n") if ref.strip()]
    elif ', ' in refs:
        ref_list = [ref.strip() for ref in refs.split(",") if ref.strip()]
    # Filter out blacklisted references
    cleaned_refs = [ref for ref in ref_list if ref.strip() not in blacklisted_refs_set]
    # Join the cleaned references back into a string

    if '\n' in refs:
        return "\n".join(cleaned_refs)
    elif ', ' in refs:
        return ", ".join(cleaned_refs)
    else:
        return "\n".join(cleaned_refs)


def get_paper_references(references):
    if not isinstance(references, str) or pd.isnull(references):
        return "-"
    # Split the 'REFERENCES' column by ';' and strip extra spaces
    reference_list = [ref.strip() for ref in references.split('\n')]
    # Filter references that are in the 'fbrf_paper_set'
    paper_references = [ref for ref in reference_list if ref in fbrf_paper_set]
    # Return the subset of references or "FALSE" if none are of type 'paper'
    return ', '.join(paper_references) if paper_references else "-"

# Defines a safe method to split column by the first space (' ')
def extract_stock_id_and_genotype(stock_value):
    if pd.isnull(stock_value):
        return pd.Series([None, None])
    stock_value = stock_value.strip()  # Remove leading/trailing spaces
    first_space_idx = stock_value.find(' ')
    if first_space_idx == -1:
        return pd.Series([stock_value, ''])
    else:
        stock_id = stock_value[:first_space_idx]
        genotype = stock_value[first_space_idx + 1:].strip()
        return pd.Series([stock_id, genotype])

def clean_genotype(genotype):
    if isinstance(genotype, str):
        return re.sub(r'\{[^}]*=([^}]*)\}', r'{\1}', genotype)
    else:
        return str(genotype)  # Return the value unchanged if it's not a string

symbol_to_name = {
    'α': 'alpha',
    'β': 'beta',
    'γ': 'gamma',
    'δ': 'delta',  # Lowercase delta
    'ε': 'epsilon',
    'ζ': 'zeta',
    'η': 'eta',
    'θ': 'theta',
    'ι': 'iota',
    'κ': 'kappa',
    'λ': 'lambda',
    'μ': 'mu',
    'ν': 'nu',
    'ξ': 'xi',
    'ο': 'omicron',
    'π': 'pi',
    'ρ': 'rho',
    'σ': 'sigma',
    'τ': 'tau',
    'υ': 'upsilon',
    'φ': 'phi',
    'χ': 'chi',
    'ψ': 'psi',
    'ω': 'omega',
    'Α': 'Alpha',  # Uppercase Alpha
    'Β': 'Beta',
    'Γ': 'Gamma',
    'Δ': 'Delta',  # Uppercase Delta
    'Ε': 'Epsilon',
    'Ζ': 'Zeta',
    'Η': 'Eta',
    'Θ': 'Theta',
    'Ι': 'Iota',
    'Κ': 'Kappa',
    'Λ': 'Lambda',
    'Μ': 'Mu',
    'Ν': 'Nu',
    'Ξ': 'Xi',
    'Ο': 'Omicron',
    'Π': 'Pi',
    'Ρ': 'Rho',
    'Σ': 'Sigma',
    'Τ': 'Tau',
    'Υ': 'Upsilon',
    'Φ': 'Phi',
    'Χ': 'Chi',
    'Ψ': 'Psi',
    'Ω': 'Omega'
}

if __name__ == "__main__":
    # Read all FlyBase published precompiled tsv files


    # Read the file (assuming it's in a tab-separated format)
    directory_path = sys.argv[1]
    input_folder_name = os.path.basename(os.path.normpath(directory_path))
    parent_dir = os.path.dirname(directory_path)

    # Define the save path in the parent directory
    save_path = os.path.join(parent_dir, "Processed_Stock_IDs")
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    txt_files = glob.glob(os.path.join(directory_path, '**.txt'), recursive=True)
    print(f'{len(txt_files)} Files found within directory: {directory_path}')
    all_txt_dfs = []
    # Define custom logic for processing <newline> as actual newlines within fields
    counter = 0
    for txt_file in txt_files:
        # Convert the cleaned text file into a DataFrame
        df = pd.read_csv(txt_file, delimiter='\t', engine='python', keep_default_na=False, quoting=csv.QUOTE_NONE,encoding='utf-8')
        df['STOCKS'] = df['STOCKS'].str.split('<newline>')
        df = df.explode('STOCKS')
        df = df.applymap(lambda x: x.replace('<newline>', '\n') if isinstance(x, str) else x)

        file_name = os.path.basename(txt_file).replace('.txt', '.csv')

        # # Construct the full save path for the output file
        file_name = file_name.replace("_flybase_gene_ids_", "_")
        file_name = file_name.replace("_Allele_Stock", "_Mutant_Stocks")
        save_file_path = os.path.join(save_path, file_name)
        df[['StockID', 'Genotype']] = df['STOCKS'].apply(extract_stock_id_and_genotype)
        df.drop(columns=["STOCKS"], inplace=True)
        mutants_no_stock = df[df["StockID"] == "-"]
        df = df[df["StockID"] != "-"]
        df = df.rename(columns={'#SUBMITTED ID': 'flybase_allele_id'})
        # Replace symbols with names in the Genotype column
        for symbol, name in symbol_to_name.items():
            df['Genotype'] = df['Genotype'].str.replace(symbol, name, regex=False)

        # Remove Annotated/ Commented parts of the flybase genotype before matching/ validating stocks
        flybase_dict["stocks_FB2024_05"]["FB_Cleaned_Genotype"] = flybase_dict["stocks_FB2024_05"]["FB_genotype"].apply(clean_genotype).str.replace("! see comment", "", regex=False).str.replace(" ", "", regex=False).str.strip()
        print(df[df['StockID'] == "32449"])
        # # Combine relevant data from tsv's
        # Get the correct stock data from flybase tsv
        df["StockID+_+Genotype"] = df["StockID"].str.replace(" ", "", regex=False) + '_+_' + df["Genotype"].str.replace(" ", "", regex=False)
        flybase_dict["stocks_FB2024_05"]['stock_number+_+FB_genotype'] = flybase_dict["stocks_FB2024_05"]['stock_number'].str.replace(" ", "", regex=False) + '_+_' + flybase_dict["stocks_FB2024_05"]['FB_Cleaned_Genotype']
        df = pd.merge(df, flybase_dict["stocks_FB2024_05"][["stock_number+_+FB_genotype", "collection_short_name", "stock_type_cv", "FB_genotype", "FB_Cleaned_Genotype", "FBst"]], how= "left", left_on = "StockID+_+Genotype", right_on = "stock_number+_+FB_genotype")
        print(save_file_path)
        print(f'Number of total rows from tsv: {df.shape[0]}')
        unmatched_mutant_stocks = df[df['stock_number+_+FB_genotype'].isna()]
        df = df[df['stock_number+_+FB_genotype'].notna()] # Remove stocks where number and genotype dont align
        df.drop(columns=["stock_number+_+FB_genotype", "StockID+_+Genotype"])
        print("After removing non matching stocks:")
        print(f'Number of total rows from tsv: {df.shape[0]}\n')
        # Add the allele information
        df = pd.merge(df, flybase_dict["fbal_to_fbgn_fb_2024_05"][["AlleleID", "AlleleSymbol"]], how = "left", left_on = "flybase_allele_id", right_on = "AlleleID")
        cols = df.columns.tolist()
        cols.insert(1, cols.pop(cols.index('AlleleSymbol')))
        # Step 1: Move the specified columns to start at index

        # Define the columns to move
        cols_to_move = ['StockID', 'Genotype', 'collection_short_name', 'stock_type_cv']

        # Remove the columns to move from their current positions
        for col in cols_to_move:
            cols.remove(col)

        # Insert these columns starting at index 5 (i.e., before the 6th column, index starts at 0)
        for i, col in enumerate(cols_to_move):
            cols.insert(5 + i, col)

        # Reorder the DataFrame columns
        df = df[cols]

        # Step 2: Group by 'ASSOCIATED_GENE', and ensure rows with StockID == '-' are at the bottom of each group

        # First, define a sorting function that moves StockID == '-' to the bottom within each group
        df_non_dash = df[df['StockID'] != '-']
        df_dash = df[df['StockID'] == '-']

        # Step 3: Group both sections by 'ASSOCIATED_GENE'

        # Define a helper function to group and sort by 'ASSOCIATED_GENE'
        def group_by_gene(grouped_df):
            return grouped_df.groupby('ASSOCIATED_GENE', group_keys=False).apply(lambda x: x.sort_values(by='StockID'))

        # Apply grouping for both sections
        df_non_dash_grouped = group_by_gene(df_non_dash)
        df_dash_grouped = group_by_gene(df_dash)

        # Step 4: Concatenate both sections: non '-' group first, then '-' group at the bottom
        df = pd.concat([df_non_dash_grouped, df_dash_grouped], ignore_index=True)
        df = df.rename(columns={"PUBLICATIONS": "REFERENCES"})
        df["REFERENCES"] = df["REFERENCES"].apply(clean_references) # Remove all refs with > 50 genes
        df['ALLELE_PAPER_REFS'] = df['REFERENCES'].apply(get_paper_references)
        df.drop(columns=["CYTOLOGY", "GENBANK_DATA", "NAME", "SYMBOL", "TAGGED_WITH", "TAGS", "STOCK_NOTES_ON_AVAILABILITY", "StockID+_+Genotype", "stock_number+_+FB_genotype"], inplace=True)

        # Save the DataFrame to the correct directory
        counter+=1
        print(f'Saving file {counter}: {file_name}')
        print(f'Size of file: {df.shape[0]} \n')
        all_txt_dfs.append(df)
    # Aggregate all txt files within a directory and output a single file
    if all_txt_dfs:
        final_df = pd.concat(all_txt_dfs, ignore_index=True)

        # Save a single output file
        combined_file_path = os.path.join(save_path, f'{input_folder_name}_Stocks.csv')
        final_df.to_csv(combined_file_path, index=False)
        print(f'Saved combined file with {final_df.shape[0]} rows at: {combined_file_path}')
    else:
        print("No data found in txt files.")
