import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np
from adjustText import adjust_text

###############################################################################
# HELPER FUNCTIONS
###############################################################################

def generate_summary_for_file(processed_data, phenotype, file_name, replicated_genes):
    """
    Generate a summary of genes (that are replicated) for the given file & phenotype.
    Previously, this function filtered rows by the difference value. We have
    removed that logic since we're now plotting raw values, not differences.
    """
    # Ensure 'ext_gene' is present.
    if 'ext_gene' not in processed_data.columns:
        processed_data = processed_data.copy()
        processed_data['ext_gene'] = ""

    # Only include replicated genes for this phenotype.
    summary = processed_data[processed_data['flybase_gene_id'].isin(replicated_genes)].copy()

    # Keep only the relevant columns.
    summary = summary[['flybase_gene_id', 'ext_gene', phenotype]].copy()

    # Add the file name.
    summary['file'] = file_name
    return summary[['file', 'flybase_gene_id', 'ext_gene', phenotype]]


def dedup_agg_data(data, phenotype, replicated_hits_for_phenotype):
    """
    Deduplicate rows for each gene. Our strategy is:
       - For any gene that is replicated, we keep the instance with the *largest absolute value*,
       - For non-replicated genes, we keep the instance with the *smallest absolute value*,
         then drop duplicates so each gene only appears once in the final dataset.
    This logic matches the original approach but uses the raw phenotype values.
    """
    data['is_replicated'] = data['flybase_gene_id'].isin(replicated_hits_for_phenotype)

    replicated_data = data[data['is_replicated']].copy()
    if not replicated_data.empty:
        replicated_data = replicated_data.sort_values(
            by=phenotype, key=lambda x: x.abs(), ascending=False
        ).drop_duplicates(subset='flybase_gene_id', keep='first')

    nonrep_data = data[~data['is_replicated']].copy()
    if not nonrep_data.empty:
        nonrep_data = nonrep_data.sort_values(
            by=phenotype, key=lambda x: x.abs(), ascending=True
        ).drop_duplicates(subset='flybase_gene_id', keep='first')

    deduped_data = pd.concat([replicated_data, nonrep_data]).sort_index()
    deduped_data.drop(columns=['is_replicated'], inplace=True)
    return deduped_data

GLOBAL_CONTROLS = {
    'TRiP attP2 ctrl': {
        'Total_sleep': 957.4139400183151,
        'Latency': 30.24826152995815,
        'P-S': 11.001618734031428,
        'Period(P-S>0)': 23.621904761904755,
    },
    'TRiP attP40 ctrl': {
        'Total_sleep': 1002.971530796179,
        'Latency': 34.902469259110404,
        'P-S': 21.075515056388454,
        'Period(P-S>0)': 23.201250000000005,
    },
    'PDF>ISO31':{
        'Total_sleep': 908.64,
        'Latency': 63.6875,
        'P-S': 99.321,
        'Period(P-S>0)': 24.4642
    }
}


def generate_rank_plot(data, phenotype, replicated_hits_for_phenotype, output_path,
                       is_agg_plot=False, control_title=""):
    """
    Generate a rank plot of the raw phenotype values (no difference from control).
    The x-axis is sorted (ranked) by phenotype, and the 3 global controls are included
    in that sorting. Labels for each control are shown above the control point with
    an upward arrow.

    Only one legend entry is kept: "Screened Genes".
    """

    # --- 1) Basic safety check ---
    if phenotype not in data.columns:
        print(f'Phenotype "{phenotype}" does not exist in the data. Skipping plot.')
        return

    # Create figure
    plt.figure(figsize=(10, 6))

    # --- 2) Deduplicate / aggregator logic ---
    if is_agg_plot:
        data['abs_value'] = data[phenotype].abs()
        data['is_replicated'] = data['flybase_gene_id'].isin(replicated_hits_for_phenotype)

        def choose_row(group):
            valid = group['abs_value'].notna()
            if valid.any():
                if group['is_replicated'].any():
                    # keep the row with largest absolute value among replicated
                    replicated_group = group[(group['is_replicated']) & valid]
                    if not replicated_group.empty:
                        return replicated_group.loc[replicated_group['abs_value'].idxmax()]
                # otherwise keep row with smallest absolute value among non-replicated
                nonrep_group = group[(~group['is_replicated']) & valid]
                if not nonrep_group.empty:
                    return nonrep_group.loc[nonrep_group['abs_value'].idxmin()]
                else:
                    return group.iloc[0]
            else:
                return group.iloc[0]

        results = []
        for gene, group in data.groupby('flybase_gene_id'):
            results.append(choose_row(group))
        data_grouped = pd.DataFrame(results)
        data_grouped.drop(columns=['abs_value', 'is_replicated'], inplace=True, errors='ignore')
    else:
        # For single-file plot, just take the average if multiple rows
        numeric_cols = data.select_dtypes(include=['number']).columns
        data_grouped = data.groupby('flybase_gene_id')[numeric_cols].mean().reset_index()

    # --- 3) Insert Global Controls into data_grouped ---
    control_rows = []
    for ctrl_name, ctrl_vals in GLOBAL_CONTROLS.items():
        if phenotype in ctrl_vals:
            ctrl_value = ctrl_vals[phenotype]
            control_rows.append({
                'flybase_gene_id': ctrl_name,
                phenotype: ctrl_value
            })
    if control_rows:
        control_df = pd.DataFrame(control_rows)
        # Concatenate controls with normal data
        data_with_ctrls = pd.concat([data_grouped, control_df], ignore_index=True)
    else:
        data_with_ctrls = data_grouped

    # --- 4) Sort combined data by this phenotype (ascending) ---
    data_sorted = data_with_ctrls.sort_values(by=phenotype).reset_index(drop=True)

    # Identify which rows are from the controls
    control_names = set(GLOBAL_CONTROLS.keys())
    is_control = data_sorted['flybase_gene_id'].isin(control_names)

    # Separate normal vs. replicated vs. controls
    normal_data = data_sorted[~is_control]
    replicated_data = normal_data[normal_data['flybase_gene_id'].isin(replicated_hits_for_phenotype)]
    control_data = data_sorted[is_control]

    # --- 5) Plot normal (screened) genes ---
    # We'll add a single legend entry for these only:
    plt.scatter(
        normal_data.index,
        normal_data[phenotype],
        label='Screened Genes',  # for the legend
        color='orange'
    )

    # --- 6) Plot replicated hits (no legend label) ---
    if not replicated_data.empty:
        plt.scatter(
            replicated_data.index,
            replicated_data[phenotype],
            color='red',
            label='_nolegend_',  # this hides them from the legend
            zorder=3
        )

    # --- 7) Plot the controls (no legend label) ---
    # plt.scatter(
    #     control_data.index,
    #     control_data[phenotype],
    #     color='red',
    #     marker='X',
    #     s=100,
    #     zorder=5,
    #     label='_nolegend_'
    # )
    plt.scatter(
        control_data.index,
        control_data[phenotype],
        color='orange',  # same color as the other genes
        marker='o',  # a normal circle marker
        s=50,  # smaller than 100 if you want it to match
        zorder=2,
        label='_nolegend_'  # ensures no additional legend entry
    )

    all_y = pd.concat([
        normal_data[phenotype],
        replicated_data[phenotype],
        control_data[phenotype]
    ])

    # 2) Compute an offset based on the total y-range:
    y_min, y_max = all_y.min(), all_y.max()
    data_range = y_max - y_min
    arrow_offset = 0.02 * data_range  # for example, 1% of the y-range

    # Annotate each control *above* its point with an upward arrow
    for i, row in control_data.iterrows():
        ctrl_name = row['flybase_gene_id']
        ctrl_value = row[phenotype]

        # We'll offset the label a bit above the point:
        # Use textcoords='offset points' for a small shift
        plt.annotate(
            ctrl_name,
            xy=(i, ctrl_value + arrow_offset),              # arrow tip at the control point
            xytext=(0, 20),                  # shift upwards by 10 points
            textcoords='offset points',
            ha='center',
            va='bottom',
            color='black',
            fontweight='bold',
            arrowprops=dict(
                width =2.5,
                color='red'
            )
        )

    # --- 8) Title ---
    title_text = phenotype
    if control_title:
        title_text += f' - {control_title}'
    plt.title("")

    # --- 9) X-axis label ---
    plt.xlabel('Genes', fontweight='bold')
    plt.xticks([])

    # --- 10) Y-axis label with units ---
    if phenotype == 'Total_sleep':
        plt.ylabel('Total Sleep time/ day (min)', fontweight='bold')
    elif phenotype == 'Latency':
        plt.ylabel('Latency (min)', fontweight='bold')
    elif phenotype == 'P-S':
        plt.ylabel('Rhythmicity (P-S)', fontweight='bold')
    elif phenotype == 'Period(P-S>0)':
        plt.ylabel('Period (P-S>0) (hrs)', fontweight='bold')
    else:
        plt.ylabel(phenotype, fontweight='bold')

    # Bold numbers on y-axis
    plt.yticks(fontweight='bold')

    # Only show one legend entry for "Screened Genes"
    plt.legend(frameon=False)

    # Turn off grid
    plt.grid(False)

    # Remove top/right spines; thicken bottom/left
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    all_y = pd.concat([
        normal_data[phenotype],
        replicated_data[phenotype],
        control_data[phenotype]
    ])

    current_ymin, current_ymax = all_y.min(), all_y.max()
    margin = 0.05 * (current_ymax - current_ymin)  # 5% margin

    plt.ylim(current_ymin - margin, current_ymax + margin)
    plt.savefig(output_path)
    plt.close()



###############################################################################
# MAIN SCRIPT
###############################################################################

def main(input_dir):
    """
    Main pipeline:
      1) Load replicated hits & relevant genes
      2) For each Excel file in the input_dir, read the raw phenotype values
      3) Generate rank plots for each phenotype (file-by-file)
      4) Aggregate across all files and generate aggregated rank plots
      5) Write out summary tables
    """
    output_dir = os.path.join(input_dir, 'Plots')
    os.makedirs(output_dir, exist_ok=True)

    agg_output_dir = os.path.join(output_dir, 'agg_plots')
    os.makedirs(agg_output_dir, exist_ok=True)

    # Load replicated hits file.
    replicated_hits_file = os.path.join(input_dir, 'replicated-hits_87-Genes.csv')
    replicated_hits_df = pd.read_csv(replicated_hits_file)
    if not {'flybase_gene_id', 'phenotype'}.issubset(replicated_hits_df.columns):
        print('Required columns ("flybase_gene_id", "phenotype") not found in replicated hits file.')
        return

    # Load relevant genes file.
    relevant_genes_file = os.path.join(input_dir, 'relevant_genes.csv')
    relevant_genes_df = pd.read_csv(relevant_genes_file)
    if 'flybase_gene_id' not in relevant_genes_df.columns:
        print('The relevant_genes.csv file is missing the column "flybase_gene_id".')
        return

    # Filter the hits to only relevant genes.
    relevant_gene_ids = set(relevant_genes_df['flybase_gene_id'])
    replicated_hits = {
        'Total_sleep': set(replicated_hits_df.loc[
            (replicated_hits_df['phenotype'] == 'Total Sleep') &
            (replicated_hits_df['flybase_gene_id'].isin(relevant_gene_ids)),
            'flybase_gene_id'
        ]),
        'Latency': set(replicated_hits_df.loc[
            (replicated_hits_df['phenotype'] == 'Latency') &
            (replicated_hits_df['flybase_gene_id'].isin(relevant_gene_ids)),
            'flybase_gene_id'
        ]),
        'P-S': set(replicated_hits_df.loc[
            (replicated_hits_df['phenotype'] == 'Rhythmicity') &
            (replicated_hits_df['flybase_gene_id'].isin(relevant_gene_ids)),
            'flybase_gene_id'
        ]),
        'Period(P-S>0)': set(replicated_hits_df.loc[
            (replicated_hits_df['phenotype'] == 'Period') &
            (replicated_hits_df['flybase_gene_id'].isin(relevant_gene_ids)),
            'flybase_gene_id'
        ])
    }

    # Define the 4 phenotype column names we expect in the Excel files:
    phenotypes = ['Total_sleep', 'Latency', 'P-S', 'Period(P-S>0)']

    # Prepare aggregated data storage
    agg_data = {ph: pd.DataFrame() for ph in phenotypes}
    summary_tables = {ph: [] for ph in phenotypes}

    # Locate Excel files
    excel_files = glob.glob(os.path.join(input_dir, '*.xlsx'))

    # Read the files and make per-file plots
    for file_path in excel_files:
        df = pd.read_excel(file_path, engine='openpyxl')

        # "Screened rows" = those with a valid FlyBase ID (FBgn)
        processed_data = df[df['flybase_gene_id'].str.startswith('FBgn', na=False)].copy()
        if processed_data.empty:
            continue

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        print(f'Processing file: {file_name}')

        # Create a directory for the plots of this file
        file_output_dir = os.path.join(output_dir, file_name)
        os.makedirs(file_output_dir, exist_ok=True)

        # For each phenotype, generate the summary and plot
        for ph in phenotypes:
            if ph not in processed_data.columns:
                continue

            # Summaries for each file
            summary = generate_summary_for_file(
                processed_data, ph, file_name, replicated_hits.get(ph, set())
            )
            if not summary.empty:
                summary_tables[ph].append(summary)

            # Generate the rank plot for each phenotype (raw values, no difference)
            # We pass the data as-is to generate_rank_plot
            plot_file_path = os.path.join(file_output_dir, f'{ph}_Rank_Plot.png')
            generate_rank_plot(processed_data, ph, replicated_hits.get(ph, set()), plot_file_path)

            # Also keep track of the data for aggregated plot
            sub = processed_data[['flybase_gene_id', ph]].copy()
            # If needed, preserve other info (like control_type) or external gene name, etc.
            # For now, we just keep gene ID + phenotype value.
            agg_data[ph] = pd.concat([agg_data[ph], sub])

    # Build aggregated plots across all files
    for ph, data in agg_data.items():
        if data.empty:
            continue

        # Deduplicate according to replicated vs non-replicated logic
        deduped_data = dedup_agg_data(data, ph, replicated_hits.get(ph, set()))

        # Output aggregated plot
        plot_file_path = os.path.join(agg_output_dir, f'{ph}_Aggregated_Rank_Plot.png')
        generate_rank_plot(deduped_data, ph, replicated_hits.get(ph, set()), plot_file_path,
                           is_agg_plot=True)

    # Save summary CSV files for each phenotype
    for ph, summaries in summary_tables.items():
        if summaries:
            combined_summary = pd.concat(summaries, ignore_index=True)
            output_csv = os.path.join(agg_output_dir, f'{ph}_summary.csv')
            combined_summary.to_csv(output_csv, index=False)
            print(f"Saved summary for {ph} to {output_csv}")

    print('All plots saved successfully!')


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Usage: python3 Phenotype_Rank_Plot_Script.py <InputDirectory>')
    else:
        main(sys.argv[1])
