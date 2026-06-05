#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

# ===============================================================
# CONFIG — edit these paths for your environment.
# ===============================================================
OUTPUT_DIR = "results/exact_exon_analysis"

# Master per-transcript processed CSV (produced by
# final_process_all_isoseq_data.py).
processed_data_path = os.path.join(OUTPUT_DIR, "all_transcripts_processed.csv")

# Where the plots are written.
base_output_dir = os.path.join(OUTPUT_DIR, "tissue_plots")
# ===============================================================
os.makedirs(base_output_dir, exist_ok=True)

# Directory specifically for sample variability plots
sample_variability_dir = os.path.join(base_output_dir, "sample_variability_plots")
os.makedirs(sample_variability_dir, exist_ok=True)

def plot_sample_variability():
    """
    Create bar plots showing isoform distribution across samples, grouped by ancestry-haplotype.
    """
    print("\nCreating sample variability bar plots...")
    
    # Load the data
    print(f"Loading data from {processed_data_path}...")
    df = pd.read_csv(processed_data_path)
    
    # Basic info
    print(f"Loaded data with {len(df)} rows and {len(df.columns)} columns")
    
    # The processed data is wild-type (WT) only.
    wt_df = df.copy()

    # Define the order for ancestry-haplotype groups
    ancestry_haplotype_order = ["EUR_H1H1", "EUR_H1H2", "AA_H1H1", "AA_H1H1_A", "AA_H1H2"]
    
    # Define the specific tissues to plot
    tissues_to_plot = {
        'Brain_MFG': 'Brain MFG',
        'Brain_Caudate': 'Brain Caudate',
        'Forebrain_neurons': 'Neurons',
        'Astrocytes': 'Astrocytes'
    }
    
    # Define isoforms of interest
    isoforms_of_interest = ["0N34_4aL", "0N3R_6", "0N4R_4aL", "0N4R_6", "1N3R", "1N3R_4aL",
                            "1N3R_4aL_6", "1N3R_6", "1N4R", "1N4R_4aL", "1N4R_6", "2N3R",
                            "2N3R_4aL", "2N3R_4aL_6", "2N3R_4a_6", "2N3R_6", "2N4R", "2N4R_4aL"]
    
    # Create distinct colors for each isoform using a combination of colormaps
    # Collect colors from several distinct colormaps to ensure variety
    distinct_colors = []
    
    # Add colors from various colormaps
    distinct_colors.extend(plt.cm.tab20.colors)  # 20 distinct colors
    distinct_colors.extend(plt.cm.tab20b.colors)  # Another 20 distinct colors
    distinct_colors.extend(plt.cm.tab20c.colors)  # Another 20 distinct colors
    
    # Add some more vibrant colors if needed
    more_colors = [
        '#FF1493',  # deep pink
        '#00FFFF',  # cyan
        '#FF4500',  # orange red
        '#9400D3',  # dark violet
        '#00FF00',  # lime
        '#1E90FF',  # dodger blue
        '#FFD700',  # gold
        '#FF00FF',  # magenta
        '#7CFC00',  # lawn green
        '#00BFFF',  # deep sky blue
    ]
    distinct_colors.extend(more_colors)
    
    # Ensure we have enough distinct colors for all isoforms
    assert len(distinct_colors) >= len(isoforms_of_interest), "Not enough distinct colors for all isoforms"
    
    # Assign a unique color to each isoform
    isoform_colors = {isoform: distinct_colors[i] for i, isoform in enumerate(isoforms_of_interest)}
    
    # Create plots for each tissue
    for tissue_code, tissue_name in tissues_to_plot.items():
        print(f"\nProcessing {tissue_name} (code: {tissue_code})...")
        
        # Check which column contains this tissue name
        tissue_column = None
        if 'Cell_Type' in wt_df.columns and tissue_code in wt_df['Cell_Type'].values:
            tissue_column = 'Cell_Type'
            tissue_df = wt_df[wt_df['Cell_Type'] == tissue_code]
        elif 'Derived_Cell_Type' in wt_df.columns and tissue_code in wt_df['Derived_Cell_Type'].values:
            tissue_column = 'Derived_Cell_Type'
            tissue_df = wt_df[wt_df['Derived_Cell_Type'] == tissue_code]
        else:
            print(f"No data found for {tissue_code}, skipping")
            continue
            
        print(f"Using {tissue_column} column for tissue type")
        
        # Check if we have enough data
        if len(tissue_df) == 0:
            print(f"No data for {tissue_name}, skipping")
            continue
        
        # Get list of samples in this tissue
        unique_samples = tissue_df['sample'].unique()
        print(f"Found {len(unique_samples)} unique samples in {tissue_name}")
        
        # Count samples in each ancestry-haplotype group to determine figure width
        group_sample_counts = {}
        for group in ancestry_haplotype_order:
            group_samples = tissue_df[tissue_df['ancestry_haplotype'] == group]['sample'].unique()
            group_sample_counts[group] = len(group_samples)
            print(f"  - {group}: {len(group_samples)} samples")
        
        # Calculate figure width based on number of samples
        # Set min width per sample and gap between groups
        min_width_per_sample = 0.6
        group_gap = 2
        
        # Calculate total width needed
        total_samples = sum(group_sample_counts.values())
        total_groups = sum(1 for count in group_sample_counts.values() if count > 0)
        
        # Calculate figure width: (samples * width_per_sample) + (gaps between groups)
        fig_width = (total_samples * min_width_per_sample) + (total_groups - 1) * 0.5 * group_gap
        
        # Ensure figure isn't too wide but is at least 10 inches
        # Make it wider to accommodate the right-side legend
        fig_width = max(min(fig_width, 32), 10)  # Between 10 and 32 inches
        
        print(f"Creating figure with width {fig_width} inches")
        
        # Create a figure with calculated size - add extra width for legend
        plt.figure(figsize=(fig_width + 5, 8))  # Add 5 inches for legend
        
        # Track positions for bars and legend
        overall_x_position = 0
        x_positions = []
        sample_ticks = []
        group_labels = []
        group_positions = []
        
        # Process each ancestry-haplotype group
        for group in ancestry_haplotype_order:
            # Filter for this group
            group_df = tissue_df[tissue_df['ancestry_haplotype'] == group]
            
            if len(group_df) == 0:
                continue
                
            # Get unique samples in this group
            samples = sorted(group_df['sample'].unique())
            
            if len(samples) == 0:
                continue
                
            # Record the center position for this group's label
            group_positions.append(overall_x_position + len(samples)/2)
            group_labels.append(group)
            
            # Process each sample
            for sample_idx, sample in enumerate(samples):
                sample_df = group_df[group_df['sample'] == sample]
                
                # Record x position for this sample
                x_positions.append(overall_x_position)
                sample_ticks.append(sample)

                # Get the total read count for this sample
                total_read_count = sample_df['total_read_count_sample'].iloc[0] if not sample_df.empty else 0

                # Track the bottom position for stacking
                bottom = 0

                # Draw stacked bars for each isoform present in this sample
                for isoform in isoforms_of_interest:
                    isoform_data = sample_df[sample_df['Isoform'] == isoform]
                    
                    if len(isoform_data) > 0:
                        # Use actual Read_Count value
                        read_count = isoform_data['Read_Count'].sum()
                        
                        # Normalize by total read count to get the proportion
                        proportion = read_count / total_read_count if total_read_count > 0 else 0
                        
                        # No need to add to legend here, we'll create a custom legend later
                        plt.bar(overall_x_position, proportion, width=0.8, 
                                color=isoform_colors[isoform], alpha=0.9, bottom=bottom)
                        
                        # Update the bottom for the next stack segment
                        bottom += proportion
                
                # Move to the next position
                overall_x_position += 1
            
            # Add a vertical line to separate groups (except after the last group)
            if group != ancestry_haplotype_order[-1] and len(samples) > 0:
                plt.axvline(x=overall_x_position - 0.5 + group_gap/2, color='black', linestyle='--', alpha=0.3)
            
            # Add some space between groups
            overall_x_position += group_gap
        
        # Remove x-axis ticks (too many samples to show labels clearly)
        plt.xticks([])
        
        # Add group labels below the x-axis
        for pos, label in zip(group_positions, group_labels):
            plt.text(pos, -0.05, label, ha='center', va='top', transform=plt.gca().get_xaxis_transform(),
                    fontsize=12, fontweight='bold')
        
        # Set y-axis limit to 1.0
        plt.ylim(0, 1.0)
        
        # Add title and labels
        plt.title(f'Isoform Distribution by Sample - {tissue_name}', fontsize=16)
        plt.ylabel('Proportion of Total Reads', fontsize=14)
        plt.xlabel('Samples by Ancestry-Haplotype Group', fontsize=14)
        
        # Add grid lines for readability
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Create a custom legend on the right side of the plot
        # Organize isoforms into two columns for better space usage
        handles = []
        
        # Add patch for each isoform
        for isoform in isoforms_of_interest:
            handles.append(mpatches.Patch(color=isoform_colors[isoform], label=isoform))
        
        # Calculate number of columns based on number of isoforms
        ncol = 1 if len(isoforms_of_interest) <= 10 else 2
        
        # Place legend on the right side
        plt.legend(handles=handles, loc='center right', bbox_to_anchor=(1.15, 0.5), 
                  fontsize=9, frameon=True, framealpha=0.8, ncol=ncol)
        
        # Add more space for the legend to prevent overlap with the figure
        plt.subplots_adjust(left=0.1, right=0.88, top=0.95, bottom=0.15)
        
        # Save the plot with a slightly lower DPI to keep file size reasonable
        plot_path = os.path.join(sample_variability_dir, f"sample_variability_{tissue_code}.png")
        
        # Check if figure dimensions are valid before saving
        fig = plt.gcf()
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        dpi = fig.dpi if hasattr(fig, 'dpi') else 100
        width_pixels = fig_width_inches * dpi
        height_pixels = fig_height_inches * dpi
        
        if width_pixels > 65000 or height_pixels > 65000:
            print(f"Warning: Figure size ({width_pixels}x{height_pixels} pixels) exceeds maximum allowed dimensions.")
            print(f"Reducing DPI to ensure figure can be saved.")
            # Calculate safe DPI
            safe_dpi = min(int(65000 / max(fig_width_inches, 1)), int(65000 / max(fig_height_inches, 1)))
            safe_dpi = max(72, min(safe_dpi, 150))  # Keep between 72 and 150 DPI
            plt.savefig(plot_path, dpi=safe_dpi, bbox_inches='tight')
        else:
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            
        plt.close()
        print(f"Saved plot to {plot_path}")

if __name__ == "__main__":
    print(f"Starting sample variability analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    plot_sample_variability()
    print(f"\nAnalysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"All outputs saved to {sample_variability_dir}")