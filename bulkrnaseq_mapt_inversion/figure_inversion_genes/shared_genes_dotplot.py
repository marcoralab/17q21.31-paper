import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ===============================================================
# CONFIG — edit for your environment.
# ===============================================================
# Directory holding the inversion-locus-filtered MFG comparison CSVs
# (output of filter_to_inversion_genes.R, with the 2 added MFG samples).
base_dir = "path/to/inversion_filtered_csvs"
# Where plots are written.
output_dir = "results/plots"
os.makedirs(output_dir, exist_ok=True)
# ===============================================================

# MFG comparisons (inversion-filtered CSVs from filter_to_inversion_genes.R)
files = {
    'H1H1_A AA vs H1H1 AA': f"{base_dir}/inversion_aa_h1h1a_vs_aa_h1h1_brain_mfg.csv",
    'H1H2 EUR vs H1H1 EUR': f"{base_dir}/inversion_eur_h1h2_vs_eur_h1h1_brain_mfg.csv",
    'H1H2 AA vs H1H1 AA':   f"{base_dir}/inversion_aa_h1h2_vs_aa_h1h1_brain_mfg.csv",
    'H1H1 EUR vs H1H1 AA':  f"{base_dir}/inversion_h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa.csv",
    'H1H2 EUR vs H1H2 AA':  f"{base_dir}/inversion_h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa.csv",
}

# Function to get significant genes separated by regulation direction
def get_significant_genes_by_direction(file_path, pvalue_threshold=0.05):
    df = pd.read_csv(file_path)
    sig_up = df[(df['pvalue'] < pvalue_threshold) & (df['log2FoldChange'] > 0)]['gene_name'].unique()
    sig_down = df[(df['pvalue'] < pvalue_threshold) & (df['log2FoldChange'] < 0)]['gene_name'].unique()
    return {'up': set(sig_up), 'down': set(sig_down)}

# Get significant genes from each file
sig_genes_by_file = {name: get_significant_genes_by_direction(path) 
                     for name, path in files.items()}

# Create separate plots for upregulated and downregulated genes
for direction in ['up', 'down']:
    # Get genes for this direction
    genes_by_file = {name: genes[direction] for name, genes in sig_genes_by_file.items()}
    
    # Get all unique genes for this direction
    all_genes = sorted(set().union(*genes_by_file.values()))
    
    # Create presence/absence matrix
    matrix = pd.DataFrame(0, index=files.keys(), columns=all_genes)
    for name, genes in genes_by_file.items():
        matrix.loc[name, list(genes)] = 1
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Create binary heatmap
    sns.heatmap(matrix, 
                cmap='Reds' if direction == 'up' else 'Blues',
                cbar=False,
                yticklabels=matrix.index,
                xticklabels=matrix.columns)
    
    # Customize the plot
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title(f'Significantly {"Upregulated" if direction == "up" else "Downregulated"} Genes Across MFG Comparisons')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_file = os.path.join(output_dir, f"shared_genes_{direction}regulated_MFG_with_flipped_geneslast2_added_2_MFG_samples.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print summary statistics
    print(f"\nSummary of {direction}regulated genes:")
    for name, genes in genes_by_file.items():
        print(f"{name}: {len(genes)} genes")
    
    print(f"\n{direction.capitalize()}regulated genes present in multiple comparisons:")
    gene_counts = matrix.sum()
    shared_genes = gene_counts[gene_counts > 1].sort_values(ascending=False)
    for gene, count in shared_genes.items():
        print(f"{gene}: present in {count} comparisons")