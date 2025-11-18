import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def main():
    # Check for command-line argument for the output directory.
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <output_directory>")
        sys.exit(1)

    output_dir = sys.argv[1]
    os.makedirs(output_dir, exist_ok=True)

    # Read in the main screen data and replicated hits file.
    df = pd.read_csv(
        "/Users/aadishms/Desktop/Projects.nosync/Allada-Lab/Neurodegeneration/Overlap_Analysis/Data/YK-Screened/RANK/genes-for-rank-plot_revised.csv")
    replicated_df = pd.read_csv(
        "/Users/aadishms/Desktop/Projects.nosync/Allada-Lab/Neurodegeneration/Overlap_Analysis/Data/YK-Screened/ManuallyProcessed/replicated-hits_87-Genes.csv")

    # Clean up extra whitespace.
    df['flybase_gene_id'] = df['flybase_gene_id'].astype(str).str.strip()
    df['ext_gene'] = df['ext_gene'].astype(str).str.strip()
    replicated_df['flybase_gene_id'] = replicated_df['flybase_gene_id'].astype(str).str.strip()
    replicated_df['phenotype'] = replicated_df['phenotype'].str.strip()

    # Mapping from phenotype key to the relevant column in df.
    phenotype_map = {
        'period': 'Period(P-S>0)',
        'latency': 'Latency',
        'total_sleep': 'Total_sleep',
        'rhythmicity': 'P-S'
    }
    # Mapping from phenotype key to the corresponding phenotype value in replicated_df.
    replicated_phenotype_map = {
        'period': 'Period',
        'latency': 'Latency',
        'total_sleep': 'Total Sleep',
        'rhythmicity': 'Rhythmicity'
    }

    # Set opacities for extreme distinction.
    alpha_control = 0.2         # Very transparent
    alpha_experimental = 0.6    # Moderately opaque
    alpha_replicated = 1.0      # Fully opaque

    # Loop over each phenotype to produce its plot.
    for phenotype_name, col_name in phenotype_map.items():
        if col_name not in df.columns:
            print(f"Column '{col_name}' not found in DataFrame. Skipping {phenotype_name}...")
            continue

        # --- Step 1: Filter & sort by the phenotype score (ascending) ---
        df_sub = df[df[col_name].notna()].copy()
        df_sub.sort_values(by=col_name, ascending=True, inplace=True)
        df_sub.reset_index(drop=True, inplace=True)
        df_sub['Rank'] = df_sub.index + 1  # x-axis: sorted order

        # --- Step 2: Identify blocks via control rows ---
        # A control row: flybase_gene_id equals '-' and ext_gene is valid (not NA and not '-').
        control_mask = (df_sub['flybase_gene_id'] == '-') & (df_sub['ext_gene'].notna()) & (df_sub['ext_gene'] != '-')
        control_indices = df_sub.index[control_mask].tolist()
        control_indices.sort()

        # Assign block numbers based on the sorted order.
        df_sub['block'] = None
        block_id = 0
        prev = 0
        for ci in control_indices:
            df_sub.loc[prev:ci, 'block'] = block_id
            block_id += 1
            prev = ci + 1
        # Keep only rows that were assigned a block.
        df_block = df_sub[df_sub['block'].notna()].copy()

        # --- Step 3: Assign group labels within each block ---
        rep_val = replicated_phenotype_map[phenotype_name]
        replicated_subset = replicated_df[replicated_df['phenotype'] == rep_val]
        replicated_ids = set(replicated_subset['flybase_gene_id'].dropna())

        def assign_group(row):
            if (row['flybase_gene_id'] == '-') and (row['ext_gene'] != '-' and pd.notna(row['ext_gene'])):
                return "control"
            elif row['flybase_gene_id'] in replicated_ids:
                return "replicated"
            else:
                return "experimental"

        df_block['group'] = df_block.apply(assign_group, axis=1)

        # --- Step 4: Plotting ---
        unique_blocks = sorted(df_block['block'].unique())
        cmap = plt.get_cmap("tab10")
        block_colors = {blk: cmap(i % 10) for i, blk in enumerate(unique_blocks)}

        fig, ax = plt.subplots(figsize=(10, 6))
        for blk in unique_blocks:
            block_df = df_block[df_block['block'] == blk]
            base_color = block_colors[blk]
            # Plot experimental genes: circles in black.
            exp_data = block_df[block_df['group'] == "experimental"]
            ax.scatter(exp_data['Rank'], exp_data[col_name], s=20, color="black",
                       alpha=alpha_experimental, marker='o')

            # Plot replicated hits: squares with black edge, in the block's base color.
            rep_data = block_df[block_df['group'] == "replicated"]
            ax.scatter(rep_data['Rank'], rep_data[col_name], s=20, color=base_color,
                       alpha=alpha_replicated, marker='s', edgecolors='black', linewidth=0.5)

            # Plot controls: triangles with black edge, in the block's base color.
            ctrl_data = block_df[block_df['group'] == "control"]
            ax.scatter(ctrl_data['Rank'], ctrl_data[col_name], s=20, color=base_color,
                       alpha=alpha_control, marker='^', edgecolors='black', linewidth=0.5)

            # Annotate controls with their ext_gene label.
            # Place the label 10 pts above and 5 pts to the left (-5 in x) of the control marker.
            for _, row in ctrl_data.iterrows():
                ax.annotate(
                    row['ext_gene'],
                    xy=(row['Rank'], row[col_name]),
                    xytext=(-10, 30),
                    textcoords="offset points",
                    fontsize=6,
                    rotation=0,
                    clip_on=True,
                    color=base_color,
                    arrowprops=dict(
                        arrowstyle="->",
                        color=base_color,
                        lw=0.5,
                        shrinkA=5,
                        shrinkB=5
                    ),
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.5)
                )

        ax.set_xlabel("Rank (sorted by " + col_name + ")")
        ax.set_ylabel(col_name)
        ax.set_title(f"Rank Plot for {phenotype_name}")

        # Create a custom legend using dummy markers.
        legend_elements = [
            Line2D([0], [0], marker='^', color='w', markerfacecolor=base_color, markersize=8,
                   markeredgecolor='black', label='Control', alpha=alpha_control),
            Line2D([0], [0], marker='o', color='w', markerfacecolor="black", markersize=8,
                   label='Experimental', alpha=alpha_experimental),
            Line2D([0], [0], marker='s', color='w', markerfacecolor=base_color, markersize=8,
                   markeredgecolor='black', label='Replicated', alpha=alpha_replicated)
        ]
        ax.legend(handles=legend_elements, title="Group")

        outpath = os.path.join(output_dir, f"{phenotype_name}_rank_plot.png")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close(fig)

    print("Plots saved successfully.")

if __name__ == "__main__":
    main()
