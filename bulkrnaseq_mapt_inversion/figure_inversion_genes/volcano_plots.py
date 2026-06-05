import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os
import logging
from datetime import datetime

# ===============================================================
# CONFIG — edit for your environment.
# ===============================================================
# Directory holding the inversion-locus-filtered MFG comparison CSVs
# (output of filter_to_inversion_genes.R, with the 2 added MFG samples).
data_dir = "path/to/inversion_filtered_csvs"
# Where volcano plots are written.
output_dir = "results/plots/volcano"
os.makedirs(output_dir, exist_ok=True)
# ===============================================================

# Set up logging
log_dir = os.path.join(output_dir, "logs")
os.makedirs(log_dir, exist_ok=True)

# Create a unique log file name with timestamp
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
log_file = os.path.join(log_dir, f'volcano_plot_inversion_{timestamp}.log')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

# The MFG comparisons plotted in the paper (inversion-filtered CSV -> plot title)
file_paths_and_titles = {
    os.path.join(data_dir, "inversion_aa_h1h1a_vs_aa_h1h1_brain_mfg.csv"): "AA H1H1A vs AA H1H1 Brain MFG",
    os.path.join(data_dir, "inversion_aa_h1h2_vs_aa_h1h1_brain_mfg.csv"): "AA H1H2 vs AA H1H1 Brain MFG",
    os.path.join(data_dir, "inversion_h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa.csv"): "H1H1 EUR vs AA Brain MFG",
    os.path.join(data_dir, "inversion_h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa.csv"): "H1H2 EUR vs AA Brain MFG",
}

# Parameters
alpha = 0.05
log2fc_threshold = 0
label_fc_threshold = 0

def create_volcano_plot(file_path, title):
    logging.info(f"Starting to process: {title}")
    
    try:
        # Read and process data
        logging.info(f"Reading file: {file_path}")
        df = pd.read_csv(file_path)
        logging.info(f"Successfully read file. Shape: {df.shape}")
        
        df['-log10(pvalue)'] = -np.log10(df['pvalue'])
        df['-log10(padj)'] = -np.log10(df['padj'])
        
        # Create plots for both p-values and adjusted p-values
        for p_type in ['pvalue', 'padj']:
            logging.info(f"Creating {p_type} plot for {title}")
            
            plt.figure(figsize=(15, 10))
            
            # Create scatter plot
            plt.scatter(df['log2FoldChange'], df[f'-log10({p_type})'],
                       c=np.where((df[p_type] < alpha) & (df['log2FoldChange'] > log2fc_threshold), 'red',
                                 np.where((df[p_type] < alpha) & (df['log2FoldChange'] < -log2fc_threshold), 'blue', 'grey')),
                       alpha=0.6, s=15)

            # Add threshold lines
            plt.axhline(-np.log10(alpha), color='black', linewidth=0.5, linestyle='--')
            plt.axvline(log2fc_threshold, color='black', linewidth=0.5, linestyle='--')
            plt.axvline(-log2fc_threshold, color='black', linewidth=0.5, linestyle='--')

            # Filter significant genes
            sig_genes = df[df[p_type] < alpha].copy()
            logging.info(f"Number of significant genes ({p_type}): {len(sig_genes)}")
            
            # Sort by significance
            sig_genes['score'] = -np.log10(sig_genes[p_type]) * abs(sig_genes['log2FoldChange'])
            sig_genes = sig_genes.sort_values('score', ascending=False)

            # Add labels
            texts = []
            for _, row in sig_genes.iterrows():
                texts.append(plt.text(row['log2FoldChange'], 
                                    row[f'-log10({p_type})'], 
                                    row['gene_name'],  # Changed from Gene_Symbol to gene_name
                                    fontsize=15,
                                    alpha=0.7))

            # Adjust text positions
            if texts:
                logging.info(f"Adjusting positions for {len(texts)} labels")
                adjust_text(texts,
                           arrowprops=dict(arrowstyle='->', color='black', lw=0.5, alpha=0.5),
                           expand_points=(1.2, 1.2),
                           force_points=(0.1, 0.1),
                           force_text=(0.2, 0.2),
                           min_arrow_dist=0.1,
                           max_time=2)

            # Customize plot
            plt.xlabel('Log2 Fold Change', fontsize=12)
            plt.ylabel(f'-Log10({p_type})', fontsize=12)
            plt.title(f'Volcano Plot: {title} (Inversion Genes)', fontsize=14, pad=20)
            plt.grid(True, alpha=0.3)

            # Add legend
            plt.legend(['Upregulated', 'Downregulated', 'Not significant'],
                      loc='upper right',
                      bbox_to_anchor=(1.15, 1))

            # Save plot
            output_file = os.path.join(output_dir, f'{title.replace(" ", "_")}_{p_type}.png')
            logging.info(f"Saving plot to: {output_file}")
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.5)
            
            # Clear memory
            plt.close('all')
            
        logging.info(f"Successfully completed processing for {title}")
        
    except Exception as e:
        logging.error(f"Error processing {title}: {str(e)}", exc_info=True)
        raise

# Start script execution
logging.info("Starting volcano plot script execution")
logging.info(f"Total number of files to process: {len(file_paths_and_titles)}")

# Process all files
successful_plots = 0
failed_plots = 0

for file_path, title in file_paths_and_titles.items():
    try:
        create_volcano_plot(file_path, title)
        successful_plots += 1
    except Exception as e:
        failed_plots += 1
        logging.error(f"Failed to process {title}", exc_info=True)

# Log summary statistics
logging.info(f"\nExecution Summary:")
logging.info(f"Total files processed: {len(file_paths_and_titles)}")
logging.info(f"Successful plots: {successful_plots}")
logging.info(f"Failed plots: {failed_plots}")
logging.info("Script execution completed")