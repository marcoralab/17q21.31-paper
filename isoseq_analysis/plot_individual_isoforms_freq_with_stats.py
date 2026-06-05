import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

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

# Define the order for ancestry-haplotype groups in plots
category_orders = {
    "WT_EUR_AA": ["EUR_H1H1", "EUR_H1H2", "AA_H1H1", "AA_H1H1_A", "AA_H1H2"],
}

# The single analysis: WT samples from EUR and AA ancestry
batches = [
    {
        "name": "WT_EUR_AA",
        "description": "WT samples from EUR and AA ancestry",
        "filters": {
            "ancestry": ["EUR", "AA"],
            "haplotype": ["H1H1", "H1H1_A", "H1H2"]
        }
    },
]

# Load the data
print(f"Loading data from {processed_data_path}...")
df = pd.read_csv(processed_data_path)

# Display basic info
print(f"Loaded data with {len(df)} rows and {len(df.columns)} columns")
print(f"Sample distribution: {df['sample'].nunique()} unique samples")
if 'Isoform' in df.columns:
    print(f"Isoform distribution: {df['Isoform'].nunique()} unique isoforms")
    # Count non-NaN isoforms
    isoform_count = df['Isoform'].notna().sum()
    print(f"Transcripts with mapped isoforms: {isoform_count} ({isoform_count/len(df)*100:.2f}%)")

# Helper function to apply batch filters
def filter_data_for_batch(df, batch_filters):
    """Apply the specified filters to the DataFrame for a batch"""
    filtered_df = df.copy()

    for column, values in batch_filters.items():
        # Include only rows whose column value is in the allowed list
        filtered_df = filtered_df[filtered_df[column].isin(values)]

    return filtered_df

# Helper function to clean up the dataframe for analysis
def prepare_data_for_analysis(df):
    # Make a copy to avoid modifying the original
    analysis_df = df.copy()
    
    # Convert frequency columns to numeric if they aren't already
    # Convert frequency columns to numeric if they aren't already
    frequency_columns = ['Raw_Frequency', 'Normalized_Frequency', 'Final_Frequency', 'read_count_over_total_sample_read_count']
    for col in frequency_columns:
        if col in analysis_df.columns:
            analysis_df[col] = pd.to_numeric(analysis_df[col], errors='coerce')
    
    # Filter out rows without isoform mapping
    if 'Isoform' in analysis_df.columns:
        analysis_df = analysis_df[analysis_df['Isoform'].notna()]
        print(f"After filtering for mapped isoforms: {len(analysis_df)} rows")
    
    # Create the combined ancestry_haplotype column if it doesn't exist
    if 'ancestry_haplotype' not in analysis_df.columns and 'ancestry' in analysis_df.columns and 'haplotype' in analysis_df.columns:
        analysis_df['ancestry_haplotype'] = analysis_df['ancestry'] + '_' + analysis_df['haplotype']

    return analysis_df

# Function to plot isoforms by cell/tissue type
def plot_isoforms_by_tissue(analysis_df, output_dir, batch_name):
    print(f"\nCreating isoform plots by tissue/brain region for batch {batch_name}...")
    
    # Identify cell/tissue type column to use
    tissue_columns = []
    if 'Derived_Cell_Type' in analysis_df.columns:
        tissue_columns.append('Derived_Cell_Type')
    if 'Cell_Type' in analysis_df.columns:
        tissue_columns.append('Cell_Type')
    
    if not tissue_columns:
        print("No tissue/cell type columns found in the data")
        return
    
    # Use the column with the most unique values
    tissue_col = max(tissue_columns, key=lambda col: analysis_df[col].nunique())
    print(f"Using {tissue_col} as tissue/brain region column")
    
    # Get list of all isoforms
    all_isoforms = sorted(analysis_df['Isoform'].unique())
    print(f"Found {len(all_isoforms)} unique isoforms")
    
    # Get list of all tissues
    all_tissues = sorted(analysis_df[tissue_col].dropna().unique())
    print(f"Found {len(all_tissues)} unique tissues/brain regions")
    
    # Get the category order for this batch
    order = category_orders.get(batch_name, None)
    
    # Create plots for both frequency types
    for freq_col in ['Raw_Frequency', 'Normalized_Frequency', 'read_count_over_total_sample_read_count']:
        print(f"\nCreating plots using {freq_col}...")
        
        # Create a subfolder for this frequency type
        freq_output_dir = os.path.join(output_dir, freq_col)
        os.makedirs(freq_output_dir, exist_ok=True)
        
        # Ancestry-haplotype plots
        if batch_name == "WT_EUR_AA":
            print(f"\nFor {batch_name}: Creating ancestry_haplotype plots with {freq_col}...")
            
            # Create ancestry_haplotype plots
            for tissue in all_tissues:
                # Filter for this tissue
                tissue_df = analysis_df[analysis_df[tissue_col] == tissue]
                
                if len(tissue_df) == 0:
                    print(f"No data for tissue: {tissue}")
                    continue
                
                if 'ancestry_haplotype' not in tissue_df.columns or tissue_df['ancestry_haplotype'].nunique() <= 1:
                    print(f"Insufficient ancestry_haplotype data for tissue: {tissue}")
                    continue
                
                # Get ancestry_haplotype groups present in this tissue
                ancestry_haplotypes = sorted(tissue_df['ancestry_haplotype'].dropna().unique())
                
                # Create directory for this tissue
                tissue_dir = os.path.join(freq_output_dir, f"{tissue.replace(' ', '_')}")
                os.makedirs(tissue_dir, exist_ok=True)
                
                # Create a plot for each isoform
                for isoform in all_isoforms:
                    # Filter for this isoform
                    isoform_df = tissue_df[tissue_df['Isoform'] == isoform]
                    
                    if len(isoform_df) == 0:
                        continue  # Skip if no data for this isoform
                    
                    # Create boxplot with individual data points
                    plt.figure(figsize=(12, 8))
                    
                    # Filter for categories that exist in this dataset
                    if order:
                        # Only include categories that exist in the data
                        available_categories = sorted(isoform_df['ancestry_haplotype'].unique())
                        plot_order = [cat for cat in order if cat in available_categories]
                        
                        # First create the boxplot with ordering
                        ax = sns.boxplot(x='ancestry_haplotype', y=freq_col, data=isoform_df, 
                                        showfliers=False, order=plot_order)
                        
                        # Then overlay with individual data points using same order
                        sns.stripplot(x='ancestry_haplotype', y=freq_col, data=isoform_df,
                                    color='black', alpha=0.5, jitter=True, order=plot_order)
                    else:
                        # If no order is defined, use default sorting
                        ax = sns.boxplot(x='ancestry_haplotype', y=freq_col, data=isoform_df, 
                                        showfliers=False)
                        
                        sns.stripplot(x='ancestry_haplotype', y=freq_col, data=isoform_df,
                                    color='black', alpha=0.5, jitter=True)
                    
                    plt.title(f'Distribution of {isoform} {freq_col} by Ancestry-Haplotype\nTissue: {tissue}')
                    plt.ylabel(f'{freq_col}')
                    plt.xlabel('Ancestry-Haplotype')
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    
                    # Save the plot
                    plot_path = os.path.join(tissue_dir, f"ancestry_haplotype_{isoform}_boxplot.png")
                    plt.savefig(plot_path, dpi=300)
                    plt.close()
                    print(f"Saved ancestry_haplotype boxplot ({freq_col}) for {isoform} in {tissue} to {plot_path}")
                
                # Create a heatmap for all isoforms and ancestry_haplotypes
                print(f"Creating ancestry_haplotype heatmap ({freq_col}) for {tissue}...")
                
                # Create pivot table: rows=isoforms, columns=ancestry_haplotype, values=frequency
                pivot_df = tissue_df.pivot_table(
                    index='Isoform',
                    columns='ancestry_haplotype',
                    values=freq_col,
                    aggfunc='mean'
                )
                
                # Reorder columns according to predefined order (if columns exist)
                if order:
                    # Only include columns that exist in the pivot table
                    existing_columns = [col for col in order if col in pivot_df.columns]
                    # If any columns in pivot_df are not in the order list, append them at the end
                    for col in pivot_df.columns:
                        if col not in existing_columns:
                            existing_columns.append(col)
                    pivot_df = pivot_df[existing_columns]
                
                # Only create heatmap if we have data
                if not pivot_df.empty:
                    plt.figure(figsize=(max(len(ancestry_haplotypes) * 1.5, 10), max(len(all_isoforms) * 0.4, 8)))
                    
                    # Create mask for NaN values
                    mask = np.isnan(pivot_df)
                    
                    # Format based on frequency type
                    fmt = '.3f' if freq_col == 'Normalized_Frequency' else '.3f'
                    
                    # Create the heatmap with mask
                    ax = sns.heatmap(pivot_df, annot=True, fmt=fmt, cmap='YlGnBu', 
                               linewidths=0.5, mask=mask, cbar_kws={'label': f'Average {freq_col}'})
                    
                    # Improve the heatmap appearance
                    plt.title(f'Average Isoform {freq_col} by Ancestry-Haplotype\nTissue: {tissue}')
                    plt.ylabel('Isoform')
                    plt.xlabel('Ancestry-Haplotype')
                    plt.xticks(rotation=45, ha='right')
                    
                    # Save the heatmap
                    heatmap_path = os.path.join(tissue_dir, f"ancestry_haplotype_all_isoforms_heatmap.png")
                    plt.tight_layout()
                    plt.savefig(heatmap_path, dpi=300)
                    plt.close()
                    print(f"Saved ancestry_haplotype heatmap ({freq_col}) for all isoforms in {tissue} to {heatmap_path}")
        
                # --- Add this snippet for subset heatmap (ancestry_haplotype branch) ---
                subset_isoforms = ["0N34_4aL", "0N3R_6", "0N4R_4aL", "0N4R_6", "1N3R", "1N3R_4aL",
                                "1N3R_4aL_6", "1N3R_6", "1N4R", "1N4R_4aL", "1N4R_6", "2N3R",
                                "2N3R_4aL", "2N3R_4aL_6", "2N3R_4a_6", "2N3R_6", "2N4R", "2N4R_4aL"]

                # Filter pivot table to only include the subset of isoforms
                subset_pivot_df = pivot_df.loc[pivot_df.index.intersection(subset_isoforms)]
                if not subset_pivot_df.empty:
                    plt.figure(figsize=(max(len(ancestry_haplotypes) * 1.5, 10), max(len(subset_pivot_df.index) * 0.4, 8)))
                    mask = np.isnan(subset_pivot_df)
                    ax = sns.heatmap(subset_pivot_df, annot=True, fmt=fmt, cmap='YlGnBu',
                                    linewidths=0.5, mask=mask, cbar_kws={'label': f'Average {freq_col}'})
                    plt.title(f'Average (Subset) Isoform {freq_col} by Ancestry-Haplotype\nTissue: {tissue}')
                    plt.ylabel('Isoform')
                    plt.xlabel('Ancestry-Haplotype')
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    subset_heatmap_path = os.path.join(tissue_dir, f"ancestry_haplotype_subset_isoforms_heatmap.png")
                    plt.savefig(subset_heatmap_path, dpi=300)
                    plt.close()
                    print(f"Saved ancestry_haplotype heatmap ({freq_col}) for subset of isoforms in {tissue} to {subset_heatmap_path}")



def analyze_statistical_significance(analysis_df, batch_name, output_dir):
    """
    Perform statistical analysis on isoform expression across groups.
    
    Parameters:
    -----------
    analysis_df : pandas DataFrame
        The filtered DataFrame for the current batch
    batch_name : str
        Name of the current batch
    output_dir : str
        Directory to save statistical results
    """
    print(f"\nPerforming statistical analysis for {batch_name}...")
    
    # Create statistical analysis output directory
    stats_output_dir = os.path.join(output_dir, "statistical_analysis")
    os.makedirs(stats_output_dir, exist_ok=True)
    
    # Determine which column to use for grouping based on batch name
    group_col = 'ancestry_haplotype'
    
    # Determine which column to use for tissue
    tissue_col = 'Derived_Cell_Type' if 'Derived_Cell_Type' in analysis_df.columns else 'Cell_Type'
    
    # Frequency columns to analyze
    freq_cols = ['Raw_Frequency', 'Normalized_Frequency', 'read_count_over_total_sample_read_count']
    freq_cols = [col for col in freq_cols if col in analysis_df.columns]
    
    # Storage for results
    all_results = []
    significant_results = []
    
    # Get list of isoforms and tissues
    all_isoforms = sorted(analysis_df['Isoform'].dropna().unique())
    all_tissues = sorted(analysis_df[tissue_col].dropna().unique())
    
    # Define reference groups for the ancestry-haplotype comparisons
    reference_groups = {
        'EUR_H1H1': ['AA_H1H1', 'EUR_H1H2'],  # Compare EUR_H1H2 to EUR_H1H1; compare AA_H1H1 to EUR_H1H1 (ancestry effect)
        'AA_H1H1': ['AA_H1H1_A', 'AA_H1H2'],  # Compare AA variants to AA_H1H1
    }
    
    # Process each isoform and tissue combination
    for isoform in all_isoforms:
        for tissue in all_tissues:
            # Filter data for current isoform and tissue
            current_df = analysis_df[(analysis_df['Isoform'] == isoform) & 
                                   (analysis_df[tissue_col] == tissue)]
            
            # Skip if not enough data
            if len(current_df) < 3 or current_df[group_col].nunique() < 2:
                continue
            
            # For each frequency column
            for freq_col in freq_cols:
                # First calculate descriptive statistics
                desc_stats = []
                group_values = {}
                
                for group in current_df[group_col].unique():
                    group_data = current_df[current_df[group_col] == group][freq_col].dropna()
                    if len(group_data) == 0:
                        continue
                    
                    # Store values for statistical tests
                    group_values[group] = group_data.values
                    
                    # Calculate descriptive statistics
                    desc_stats.append({
                        'Isoform': isoform,
                        'Tissue': tissue,
                        'Group': group,
                        'Count': len(group_data),
                        'Mean': group_data.mean(),
                        'Median': group_data.median(),
                        'StdDev': group_data.std(),
                        'SEM': group_data.std() / np.sqrt(len(group_data)),
                        'Frequency_Type': freq_col
                    })
                
                # Save descriptive statistics
                if desc_stats:
                    desc_df = pd.DataFrame(desc_stats)
                    desc_file = os.path.join(stats_output_dir, 
                                           f"{isoform}_{tissue}_{freq_col}_descriptive_stats.csv")
                    desc_df.to_csv(desc_file, index=False)
                
                # Skip if not enough groups with data
                if len(group_values) < 2:
                    continue
                
                # Perform statistical tests
                # 1. If only 2 groups, do t-test
                if len(group_values) == 2:
                    groups = list(group_values.keys())
                    t_stat, p_value = stats.ttest_ind(
                        group_values[groups[0]], 
                        group_values[groups[1]], 
                        equal_var=False
                    )
                    
                    # Format result with explicit reference structure
                    # Determine reference and comparison groups
                    if reference_groups and groups[0] in reference_groups and groups[1] in reference_groups[groups[0]]:
                        ref_group, comp_group = groups[0], groups[1]
                    elif reference_groups and groups[1] in reference_groups and groups[0] in reference_groups[groups[1]]:
                        ref_group, comp_group = groups[1], groups[0]
                        # Negate t-stat to maintain correct direction
                        t_stat = -t_stat
                    else:
                        # Default if no reference structure defined
                        ref_group, comp_group = groups[0], groups[1]
                    
                    result = {
                        'Isoform': isoform,
                        'Tissue': tissue,
                        'Frequency_Type': freq_col,
                        'Test': 't-test',
                        'Reference_Group': ref_group,
                        'Compared_Group': comp_group,
                        'Comparison': f"{comp_group} vs {ref_group}",
                        'Statistic': t_stat,
                        'P_Value': p_value,
                        'Significant': p_value < 0.05
                    }
                    
                    all_results.append(result)
                    if p_value < 0.05:
                        significant_results.append(result)
                
                # 2. If more than 2 groups, do ANOVA + Tukey HSD
                elif len(group_values) > 2:
                    # First perform ANOVA
                    anova_data = list(group_values.values())
                    f_stat, p_value = stats.f_oneway(*anova_data)
                    
                    # Store ANOVA result
                    anova_result = {
                        'Isoform': isoform,
                        'Tissue': tissue,
                        'Frequency_Type': freq_col,
                        'Test': 'ANOVA',
                        'Groups': ", ".join(group_values.keys()),
                        'F_Statistic': f_stat,
                        'P_Value': p_value,
                        'Significant': p_value < 0.05
                    }
                    
                    all_results.append(anova_result)
                    
                    # If ANOVA is significant, do post-hoc tests
                    if p_value < 0.05:
                        significant_results.append(anova_result)
                        
                        # Prepare data for Tukey's HSD
                        all_values = []
                        all_groups = []
                        
                        for group, values in group_values.items():
                            all_values.extend(values)
                            all_groups.extend([group] * len(values))
                        
                        # Perform Tukey's HSD
                        tukey = pairwise_tukeyhsd(
                            endog=all_values,
                            groups=all_groups,
                            alpha=0.05
                        )
                        
                        # Convert to DataFrame
                        tukey_df = pd.DataFrame(
                            data=tukey._results_table.data[1:],
                            columns=tukey._results_table.data[0]
                        )
                        
                        # Save full Tukey results
                        tukey_file = os.path.join(stats_output_dir, 
                                                f"{isoform}_{tissue}_{freq_col}_tukey_results.csv")
                        tukey_df.to_csv(tukey_file, index=False)
                        
                        # Process each pairwise comparison with reference structure
                        for _, row in tukey_df.iterrows():
                            group1 = row['group1']
                            group2 = row['group2']
                            
                            

                            # Determine reference group according to reference structure
                            if reference_groups:
                                if group1 in reference_groups and group2 in reference_groups[group1]:
                                    ref_group, comp_group = group1, group2
                                    mean_diff = row['meandiff']
                                elif group2 in reference_groups and group1 in reference_groups[group2]:
                                    ref_group, comp_group = group2, group1
                                    mean_diff = -row['meandiff']  # Negate to maintain correct direction
                                else:
                                    # If no reference relationship defined, use default
                                    ref_group, comp_group = 'EUR_H1H1', group2 if group2 != 'EUR_H1H1' else group1
                                    mean_diff = row['meandiff'] if ref_group == group1 else -row['meandiff']
                            else:
                                # Default to alphabetical if no reference structure
                                ref_group, comp_group = (group1, group2) if group1 < group2 else (group2, group1)
                                mean_diff = row['meandiff'] if group1 < group2 else -row['meandiff']
                            
                            # Only add significant results to summary
                            if row['reject']:
                                tukey_result = {
                                    'Isoform': isoform,
                                    'Tissue': tissue,
                                    'Frequency_Type': freq_col,
                                    'Test': 'Tukey HSD',
                                    'Reference_Group': ref_group,
                                    'Compared_Group': comp_group,
                                    'Comparison': f"{comp_group} vs {ref_group}",
                                    'Mean_Diff': mean_diff,
                                    'P_Value': row['p-adj'],
                                    'Significant': True
                                }
                                
                                significant_results.append(tukey_result)
    
    # Save results
    if all_results:
        all_df = pd.DataFrame(all_results)
        all_file = os.path.join(stats_output_dir, f"{batch_name}_all_statistical_tests.csv")
        all_df.to_csv(all_file, index=False)
        print(f"Saved all statistical results to {all_file}")
    
    if significant_results:
        sig_df = pd.DataFrame(significant_results)
        sig_file = os.path.join(stats_output_dir, f"{batch_name}_significant_results.csv")
        sig_df.to_csv(sig_file, index=False)
        print(f"Found {len(significant_results)} statistically significant results")
        print(f"Saved significant results to {sig_file}")
        
        # Check for 1N4R in Brain_MFG as in your example
        if batch_name == "WT_EUR_AA":
            mfg_1n4r = sig_df[(sig_df['Isoform'] == '1N4R') & 
                             (sig_df['Tissue'] == 'Brain_MFG')]
            
            if not mfg_1n4r.empty:
                print("\nResults for 1N4R in Brain_MFG:")
                for _, row in mfg_1n4r.iterrows():
                    if row['Test'] == 'ANOVA':
                        print(f"  ANOVA: F={row['F_Statistic']:.3f}, p={row['P_Value']:.4f}")
                    else:
                        print(f"  {row['Comparison']}: {row['Test']}, p={row['P_Value']:.4f}")
    else:
        print("No statistically significant results found")
    
    return significant_results

# Function to create enhanced plots with significance indicators
def create_significance_plots(analysis_df, significant_results, batch_name, output_dir):
    """Create plots with significance indicators for significant results."""
    if not significant_results:
        return
    
    # Convert to DataFrame if it's a list
    if isinstance(significant_results, list):
        sig_df = pd.DataFrame(significant_results)
    else:
        sig_df = significant_results
    
    # Skip if empty
    if sig_df.empty:
        return
    
    print("\nCreating enhanced plots with significance indicators...")
    
    # Create plots directory
    plots_dir = os.path.join(output_dir, "significance_plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # Determine group and tissue columns
    group_col = 'ancestry_haplotype'
    tissue_col = 'Derived_Cell_Type' if 'Derived_Cell_Type' in analysis_df.columns else 'Cell_Type'
    
    # Get unique isoform-tissue-freqcol combinations with significant results
    combinations = sig_df[['Isoform', 'Tissue', 'Frequency_Type']].drop_duplicates()
    
    # Process each combination
    for _, combo in combinations.iterrows():
        isoform = combo['Isoform']
        tissue = combo['Tissue']
        freq_col = combo['Frequency_Type']
        
        # Get ANOVA results if any
        anova_results = sig_df[(sig_df['Isoform'] == isoform) & 
                              (sig_df['Tissue'] == tissue) & 
                              (sig_df['Frequency_Type'] == freq_col) & 
                              (sig_df['Test'] == 'ANOVA')]
        
        # Get post-hoc results if any
        posthoc_results = sig_df[(sig_df['Isoform'] == isoform) & 
                               (sig_df['Tissue'] == tissue) & 
                               (sig_df['Frequency_Type'] == freq_col) & 
                               (sig_df['Test'] != 'ANOVA')]
        
        # Filter data
        current_df = analysis_df[(analysis_df['Isoform'] == isoform) & 
                               (analysis_df[tissue_col] == tissue)]
        
        if len(current_df) == 0:
            continue
        
        # Create the plot
        plt.figure(figsize=(12, 8))
        
        # Get the category order for this batch from your existing code
        order = category_orders.get(batch_name, None)
        
        # Filter for categories that exist in this dataset
        if order:
            # Only include categories that exist in the data
            available_categories = sorted(current_df[group_col].unique())
            plot_order = [cat for cat in order if cat in available_categories]
            
            # Create the boxplot with ordering
            ax = sns.boxplot(x=group_col, y=freq_col, data=current_df, 
                            showfliers=False, order=plot_order)
            
            # Add individual data points using same order
            sns.stripplot(x=group_col, y=freq_col, data=current_df,
                        color='black', alpha=0.5, jitter=True, order=plot_order)
        else:
            # If no order is defined, use default sorting
            ax = sns.boxplot(x=group_col, y=freq_col, data=current_df, 
                            showfliers=False)
            
            sns.stripplot(x=group_col, y=freq_col, data=current_df,
                        color='black', alpha=0.5, jitter=True)
        
        # Add significance markers if we have post-hoc results
        if not posthoc_results.empty:
            # Get y-axis range for positioning significance bars
            y_min, y_max = plt.ylim()
            y_range = y_max - y_min
            start_y = y_max + (y_range * 0.05)
            bar_height = y_range * 0.05
            
            # Map groups to x-coordinates
            group_positions = {}
            for i, tick in enumerate(ax.get_xticks()):
                if i < len(plot_order if order else current_df[group_col].unique()):
                    group_name = plot_order[i] if order else sorted(current_df[group_col].unique())[i]
                    group_positions[group_name] = i
            
            # Add significance bars
            for i, (_, row) in enumerate(posthoc_results.iterrows()):
                if 'Reference_Group' in row and 'Compared_Group' in row:
                    ref_group = row['Reference_Group']
                    comp_group = row['Compared_Group']
                    
                    # Skip if any group not in positions
                    if ref_group not in group_positions or comp_group not in group_positions:
                        continue
                    
                    # Get x-coordinates
                    x1 = group_positions[ref_group]
                    x2 = group_positions[comp_group]
                    
                    # Calculate y-coordinate for this bar (higher for each additional bar)
                    y = start_y + (i * bar_height * 1.5)
                    
                    # Draw the significance bar
                    plt.plot([x1, x2], [y, y], 'k-', linewidth=1.5)
                    
                    # Add significance stars
                    p_val = row['P_Value']
                    if p_val < 0.001:
                        sig_str = '***'
                    elif p_val < 0.01:
                        sig_str = '**'
                    elif p_val < 0.05:
                        sig_str = '*'
                    else:
                        sig_str = 'ns'
                    
                    plt.text((x1+x2)/2, y + (bar_height*0.5), sig_str, 
                           ha='center', va='bottom', fontsize=12)
        
        # Add title with ANOVA results if available
        if not anova_results.empty:
            f_stat = anova_results.iloc[0]['F_Statistic']
            p_val = anova_results.iloc[0]['P_Value']
            plt.title(f'{isoform} in {tissue}\nANOVA: F={f_stat:.2f}, p={p_val:.4f}')
        else:
            plt.title(f'Distribution of {isoform} {freq_col} by {group_col}\nTissue: {tissue}')
        
        plt.ylabel(f'{freq_col}')
        plt.xlabel(group_col)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        # Save the plot
        plot_path = os.path.join(plots_dir, f"{isoform}_{tissue}_{freq_col}_significance.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        
        print(f"Saved significance plot for {isoform} in {tissue} ({freq_col})")

# Function to run analysis for a batch
def run_batch_analysis(batch_name, batch_description, batch_filters):
    print(f"\n{'='*80}")
    print(f"BATCH: {batch_name}")
    print(f"DESCRIPTION: {batch_description}")
    print(f"FILTERS: {batch_filters}")
    print(f"{'='*80}")
    
    # Create output directory for this batch
    batch_output_dir = base_output_dir
    os.makedirs(batch_output_dir, exist_ok=True)
    
    # Filter data for this batch
    batch_df = filter_data_for_batch(df, batch_filters)
    
    # Check if we have data
    if len(batch_df) == 0:
        print(f"No data found for batch {batch_name} with filters {batch_filters}")
        return
    
    # Print batch data info
    print(f"Filtered data: {len(batch_df)} rows, {batch_df['sample'].nunique()} samples")
    
    if 'ancestry' in batch_df.columns:
        print("Ancestry distribution:")
        print(batch_df['ancestry'].value_counts())
    
    if 'haplotype' in batch_df.columns:
        print("Haplotype distribution:")
        print(batch_df['haplotype'].value_counts())
    
    # Check if we have cell type or derived cell type
    cell_type_cols = []
    if 'Cell_Type' in batch_df.columns:
        cell_type_cols.append('Cell_Type')
        print("Cell Type distribution:")
        print(batch_df['Cell_Type'].value_counts())
    
    if 'Derived_Cell_Type' in batch_df.columns:
        cell_type_cols.append('Derived_Cell_Type')
        print("Derived Cell Type distribution:")
        print(batch_df['Derived_Cell_Type'].value_counts())
    
    if not cell_type_cols:
        print("No cell type or tissue columns found. Cannot proceed with tissue-based analysis.")
        return
    
    # Prepare data for analysis
    analysis_df = prepare_data_for_analysis(batch_df)
    
    if 'Isoform' in analysis_df.columns and analysis_df['Isoform'].notna().sum() > 0:
        print(f"\nStarting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Output directory: {batch_output_dir}")
        
        # Plot isoforms by tissue
        plot_isoforms_by_tissue(analysis_df, batch_output_dir, batch_name)
        
        # Run statistical analysis - ADD THIS BLOCK
        significant_results = analyze_statistical_significance(analysis_df, batch_name, batch_output_dir)
        
        # Create enhanced plots with significance markers
        if significant_results:
            create_significance_plots(analysis_df, significant_results, batch_name, batch_output_dir)
        # END OF ADDED BLOCK
        
        # Save filtered data for this batch
        batch_data_path = os.path.join(batch_output_dir, f"{batch_name}_filtered_data.csv")
        analysis_df.to_csv(batch_data_path, index=False)
        print(f"Saved filtered data to {batch_data_path}")
        
        print(f"\nAnalysis for batch {batch_name} completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    else:
        print(f"Error: No valid isoform data found in the batch {batch_name}")

# Main execution
if __name__ == "__main__":
    print(f"Starting tissue-based isoform analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Run analysis for each batch
    for batch in batches:
        run_batch_analysis(batch['name'], batch['description'], batch['filters'])
    
    print(f"\nAll batch analyses completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"All outputs saved to {base_output_dir}")