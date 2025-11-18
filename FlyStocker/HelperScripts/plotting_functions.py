import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import statsmodels.api as sm
from statsmodels.graphics.mosaicplot import mosaic
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from adjustText import adjust_text
import os


def plot_gene_stock_histogram(df, csv_file, stock_sheet_folder):
     # Count the number of occurrences (rows) for each unique gene
    gene_counts = df['ASSOCIATED_GENE'].value_counts()

    # Sort the gene counts by their values (number of stocks)
    sorted_gene_counts = gene_counts.sort_values()

    # Create an index for the x-axis, starting from 1
    gene_indices = range(1, len(sorted_gene_counts) + 1)

    # Plotting the histogram with x as the index and y as the number of stocks
    plt.figure(figsize=(12, 8))  # Increase figure size for readability
    plt.bar(gene_indices, sorted_gene_counts.values, color='blue')  # Use a single color (blue)

    # Set axis labels and title
    plt.xlabel('Gene Index (sorted by number of stock, ascending)', fontsize=14)
    plt.ylabel('Number of Stocks', fontsize=14)
    plt.title(f'Gene vs. Number of Stocks in {os.path.splitext(os.path.basename(csv_file))[0]}', fontsize=16)

    # Add gridlines for both x and y axes
    # Add gridlines for both x and y axes with finer gradations
    plt.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)  # Major gridlines
    plt.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)  # Minor gridlines

    # Add more y-axis gradations
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=15))  # Control the number of y-axis ticks
    plt.gca().yaxis.set_minor_locator(MaxNLocator(nbins=30))  # Add minor y-axis ticks for more lines

    # Reduce the number of x-tick labels to avoid overlap (e.g., show every 10th label)
    tick_step = 10  # Show every 10th label
    plt.xticks(gene_indices[::tick_step], gene_indices[::tick_step], rotation=90, fontsize=6)

    # Tighten layout to ensure the plot fits well
    plt.tight_layout()

    # Save the plot
    output_path = os.path.join(stock_sheet_folder, f'stock_histogram_{os.path.splitext(os.path.basename(csv_file))[0]}.png')
    plt.savefig(output_path, dpi=300)

def plot_cumulative_gene_histogram(df, csv_file, stock_sheet_folder, fieldOfInterest = "ASSOCIATED_GENE"):
    # Count the number of occurrences (rows) for each unique gene
    gene_counts = df[fieldOfInterest].value_counts()

    # Sort the gene counts in ascending order (so we can easily count the genes with >= x stocks)
    sorted_gene_counts = gene_counts.sort_values(ascending=True)

    # Total number of unique genes
    total_genes = len(gene_counts)

    # Create a list to store the cumulative count of genes with at least x stocks
    cumulative_genes = []
    cumulative_percentages = []

    # Iterate through possible stock counts (from 1 to max stocks)
    for threshold in range(1, sorted_gene_counts.max() + 1):
        # Count the number of genes that have at least 'threshold' number of stocks
        num_genes_at_least_x = (gene_counts >= threshold).sum()
        cumulative_genes.append(num_genes_at_least_x)

        # Calculate the percentage of total genes with at least 'threshold' stocks
        percentage = (num_genes_at_least_x / total_genes) * 100
        cumulative_percentages.append(percentage)

    # The x-axis will represent the number of stocks (threshold)
    x_values = range(1, len(cumulative_genes) + 1)

    # Create a figure and axis for the plot
    plt.figure(figsize=(12, 8))

    # Plot the cumulative gene count where y represents genes with at least x stocks
    plt.plot(x_values, cumulative_genes, marker='o', linestyle='-', color='blue', label='Genes with at least X Stocks')

    # Plot the percentages as a secondary y-axis
    ax = plt.gca()  # Get the current axis
    ax2 = ax.twinx()  # Create a second y-axis that shares the same x-axis
    ax2.plot(x_values, cumulative_percentages, marker='x', linestyle='--', color='black', label='Percentage of Genes')

    # Set axis labels and title
    ax.set_xlabel('Number of Stocks (Threshold)', fontsize=14)
    ax.set_ylabel(f'Number of {fieldOfInterest}s with at least X Stocks', fontsize=14)
    ax.set_title(f'Number of {fieldOfInterest}s with at least X Stocks for {os.path.splitext(os.path.basename(csv_file))[0]}', fontsize=16)

    # Increase the number of ticks on both x and y axes
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=15))  # More ticks on x-axis
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))  # More ticks on the primary y-axis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=10))  # More ticks on the secondary y-axis

    # Minor ticks: add extra gridlines for better granularity
    plt.minorticks_on()

    # Customize gridlines for both major and minor ticks
    ax.grid(which='major', linestyle='-', linewidth=0.75, alpha=0.8)  # Solid lines for major ticks
    ax.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.6)  # Dotted lines for minor ticks

    # Add text annotations for total number of genes and stocks
    total_stocks = gene_counts.sum()
    ax.text(0.8, 0.2, f'Total {fieldOfInterest}: {total_genes}\nTotal Stocks: {total_stocks}', fontsize=12, transform=ax.transAxes)

    # Annotate every few points on the percentage line
    for i in range(0, len(x_values), max(1, len(x_values) // 15)):  # Space out annotations
        ax2.annotate(f'{cumulative_percentages[i]:.1f}%', (x_values[i], cumulative_percentages[i]), 
                     textcoords="offset points", xytext=(5,15), ha='center', color='orange', fontsize=9)

    # Tighten layout and save the plot
    plt.tight_layout()
    histogram_type = ""
    if(fieldOfInterest == "ASSOCIATED_GENE"):
        histogram_type = "gene"
    else:
        histogram_type = fieldOfInterest
    output_path = os.path.join(stock_sheet_folder, f'stock_{histogram_type}_cumulative_distribution_{os.path.splitext(os.path.basename(csv_file))[0]}.png')
    plt.savefig(output_path, dpi=300)

def gene_map_join(df, gene_mapping_df):
    return

def insertion_join(primary_df, insertion_mapping_df, insertion_id_col='insertion_id', mapping_fbti_col='FBti#', mapping_symbol_col='insertion_symbol'):
    """
    Merges `primary_df` with `insertion_mapping_df` based on a join key.
    If `insertion_id` in `primary_df` starts with "FBti", it uses this ID; otherwise, it uses the insertion symbol.

    Parameters:
    - primary_df (pd.DataFrame): The primary DataFrame containing insertion information.
    - insertion_mapping_df (pd.DataFrame): The insertion mapping DataFrame.
    - insertion_id_col (str): Column in `primary_df` with FBti ID or insertion symbol.
    - mapping_fbti_col (str): Column in `insertion_mapping_df` with FBti IDs.
    - mapping_symbol_col (str): Column in `insertion_mapping_df` with insertion symbols.

    Returns:
    - pd.DataFrame: The resulting DataFrame after the merge.
    """
 # Step 1: Create join key in primary_df
    primary_df['join_key'] = primary_df[insertion_id_col].str.strip()

    # Step 2: Standardize insertion_mapping_df columns for matching
    insertion_mapping_df[mapping_fbti_col] = insertion_mapping_df[mapping_fbti_col].astype(str).str.strip()
    insertion_mapping_df[mapping_symbol_col] = insertion_mapping_df[mapping_symbol_col].str.strip()

    # Step 3: Split primary_df into two parts based on the join condition
    fbti_df = primary_df[primary_df[insertion_id_col].str.startswith("FBti", na=False)]
    symbol_df = primary_df[~primary_df[insertion_id_col].str.startswith("FBti", na=False)]

    # Step 4: Merge `fbti_df` with `insertion_mapping_df` on `FBti#`
    fbti_merged = fbti_df.merge(
        insertion_mapping_df,
        left_on='join_key',
        right_on=mapping_fbti_col,
        how='left'
    )

    # Step 5: Merge `symbol_df` with `insertion_mapping_df` on `insertion_symbol`
    symbol_merged = symbol_df.merge(
        insertion_mapping_df,
        left_on='join_key',
        right_on=mapping_symbol_col,
        how='left'
    )

    # Step 6: Concatenate the two merged DataFrames
    merged_df = pd.concat([fbti_merged, symbol_merged], ignore_index=True)

    # Step 7: Drop unnecessary '_x' or '_y' columns
    merged_df = merged_df.loc[:, ~merged_df.columns.str.endswith(('_x', '_y'))]
    
    return merged_df

def drop_suffix_columns(df, suffixes=('_x', '_y')):
    """
    Drops columns in the DataFrame that end with specified suffixes.

    Parameters:
    - df (pd.DataFrame): The DataFrame from which to drop columns.
    - suffixes (tuple): A tuple of suffixes to drop (default is ('_x', '_y')).

    Returns:
    - pd.DataFrame: A DataFrame with the specified columns removed.
    """
    # Select only columns that do not end with the specified suffixes
    return df.loc[:, ~df.columns.str.endswith(suffixes)]


def plot_stock_distribution(df, csv_file, stock_sheet_folder, fieldOfInterest = "ASSOCIATED_GENE"):
    if "AlleleSymbol_x" in df.columns:
        df["AlleleSymbol"] = df["AlleleSymbol_x"]
    gene_counts = df[fieldOfInterest].value_counts()

    # Count how many genes have each specific number of stocks
    stock_distribution = gene_counts.value_counts().sort_index()

    # Calculate the total number of stocks
    total_stocks = df.shape[0]
    total_genes = df[fieldOfInterest].nunique()
    # Convert the index to a numeric array (number of stocks per gene)
    stocks_per_gene = stock_distribution.index.to_numpy()
    fig, ax1 = plt.subplots(figsize=(12, 8))
    # Calculate the cumulative sum of the stocks covered by the genes at each x value
    top_10_genes = gene_counts.head(10)
    unique_colors = sns.color_palette('Set1', n_colors=10)
    bar_colors = ['gray'] * len(gene_counts)
    # Create a figure and axis for the bar plot


    top_gene_color_map = {}

    for i, (gene, count) in enumerate(top_10_genes.items()):
        bar_index = np.where(stock_distribution.index == count)[0]
        if len(bar_index) > 0:
            bar_colors[bar_index[0]] = unique_colors[i]
            top_gene_color_map[gene] = unique_colors[i]

    # Plot the bar chart for stock distribution, using the custom colors
    bars = ax1.bar(stock_distribution.index, stock_distribution.values, color=bar_colors, edgecolor=bar_colors)

    # Set labels and title with orange te
    ax1.set_xlabel(f'Number of Stocks per {fieldOfInterest}', fontsize=14, color='black')
    ax1.set_ylabel(f'Number of Unique {fieldOfInterest}', fontsize=14, color='black')
    ax1.set_title(f'Distribution of {fieldOfInterest} by Number of Stocks in {os.path.splitext(os.path.basename(csv_file))[0]}', fontsize=16, color='black')

    # Add gridlines for better readability
    ax1.grid(True, axis='y')
    ax1.grid(True, axis='x')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))  # More ticks on x-axis
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))  # More ticks on the primary y-axis

    # Minor ticks: add extra gridlines for better granularity
    plt.minorticks_on()

    # Customize gridlines for both major and minor ticks
    ax1.grid(which='major', linestyle='-', linewidth=0.75, alpha=0.8)  # Solid lines for major ticks
    ax1.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.6)  # Dotted lines for minor ticks

    # Customize the font size of the tick labels
    ax1.tick_params(axis='both', which='major', labelsize=10)  # Smaller font for major ticks
    ax1.tick_params(axis='both', which='minor', labelsize=8)   #

    # Customize x and y axis ticks, and change tick colors to orange
    ax1.tick_params(axis='x', labelsize=10, rotation=45, colors='black')
    ax1.tick_params(axis='y', labelsize=10, colors='black')
    
    # Add text annotations for total number of stocks, with orange text
    ax1.text(0.8, 0.3, f'Total Stocks: {total_stocks}', fontsize=12, transform=ax1.transAxes, color='black')
    ax1.text(0.8, 0.35, f'Total {fieldOfInterest}s: {total_genes}', fontsize=12, transform=ax1.transAxes, color='black')

   # Label the top 10 fieldOfInterest candidates on the plot as legend
    legend_labels = [f'{gene}: {count} stocks' for gene, count in top_10_genes.items()]
    legend_patches = [plt.Line2D([0], [0], marker='o', color=color, lw=0, markersize=10, label=label) 
                      for label, color in zip(legend_labels, unique_colors)]

    ax1.legend(handles=legend_patches, title=f'Top 10 {fieldOfInterest}s by Stock Count', loc='upper right', fontsize=10, title_fontsize=12)

    # Tighten layout to ensure the plot fits well
    plt.tight_layout()

    # Save the plot
    output_path = os.path.join(stock_sheet_folder, f'stock_distribution_by_{fieldOfInterest}_{os.path.splitext(os.path.basename(csv_file))[0]}.png')
    plt.savefig(output_path, dpi=300, bbox_inches = 'tight')

# TODO: Add logic to handle GENE_HAS_SLEEP_PAPER/ GENE_HAS_SLEEP_PAPER_REF col
def plot_mosaic_style_plot(stock_hits, csv_file, stock_sheet_folder):
    # Define your categories for filtering
    stock_hits['Bloomington'] = stock_hits['collection_short_name'].str.contains('Bloomington', na=False)
    stock_hits['With Paper Reference'] = stock_hits['PAPER_REFERENCES'] != '-'
    stock_hits['With UAS'] = stock_hits['NATURE_OF_LESION'].str.contains('UAS', na=False)
    stock_hits['With sgRNA'] = stock_hits['NATURE_OF_LESION'].str.contains('sgRNA', na=False)
    stock_hits['Neither UAS nor sgRNA'] = (~stock_hits['NATURE_OF_LESION'].str.contains('UAS', na=False)) & (~stock_hits['NATURE_OF_LESION'].str.contains('sgRNA', na=False))

    # Create a cross-tabulation of the categories
    crosstab = pd.crosstab(index=[stock_hits['Bloomington'], stock_hits['With Paper Reference']],
                           columns=[stock_hits['With UAS'], stock_hits['With sgRNA'], stock_hits['Neither UAS nor sgRNA']])

    # Create a figure and axis with larger size for better readability
    fig, ax = plt.subplots(figsize=(14, 10))

    # Define distinct colors for each combination of categories (mimicking a mosaic plot)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    # Plot the stacked bar chart
    crosstab.plot(kind='bar', stacked=True, ax=ax, color=colors, edgecolor='black')

    # Set axis labels and title with larger font sizes
    ax.set_title('Mosaic-style Bar Plot: UAS, sgRNA, Bloomington, and Paper Reference Categories', fontsize=16, pad=50)
    ax.set_xlabel('Bloomington Status & Paper Reference', fontsize=14)
    ax.set_ylabel('Number of Stocks', fontsize=14)

    # Get the total number of stocks for percentage calculation
    total_stocks = stock_hits.shape[0]

    # Annotate the bars with the actual counts and percentages
    # Annotate the bars with the actual counts and percentages
    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        if height > 0:
            # Label format for counts and percentages
            label = f'{height:.0f}\n({height * 100 / total_stocks:.1f}%)'
            ax.annotate(label, (x + width / 2, y + height / 2),
                        ha='center', va='center', fontsize=10, color='white', weight='bold')

    # Customize the x-tick labels to clearly represent the categories
    ax.set_xticklabels(['Bloomington & With Paper Ref.', 'Bloomington & No Paper Ref.', 
                        'Non-Bloomington & With Paper Ref.', 'Non-Bloomington & No Paper Ref.'],
                    rotation=45, ha='right', fontsize=12)

    # Add a custom legend that reflects meaningful labels instead of True/False values
    custom_legend_labels = ['With UAS', 'With sgRNA', 'No UAS, With sgRNA', 'No UAS, No sgRNA']
    legend_patches = [plt.Line2D([0], [0], color=color, lw=4) for color in colors[:len(custom_legend_labels)]]
    ax.legend(legend_patches, custom_legend_labels, title='Legend: UAS/sgRNA Categories', fontsize=12, title_fontsize=14)

    # Add text labels for each True/False combination, similar to a mosaic plot
    for i, (bloomington, paper_ref) in enumerate(crosstab.index):
        label = f'Bloomington: {"T" if bloomington else "F"}, Paper Ref: {"T" if paper_ref else "F"}'
        ax.annotate(label, xy=(i, -50), ha='center', fontsize=10, rotation=0)

    # Add a grid for clarity
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)

    # Adjust layout for better spacing
    plt.subplots_adjust(top=0.85, bottom=0.25)

    # Save the bar plot to a file
    output_path = os.path.join(stock_sheet_folder, f'stock_mosaic_style_bar_{os.path.splitext(os.path.basename(csv_file))[0]}.png')
    plt.savefig(output_path, dpi=300)