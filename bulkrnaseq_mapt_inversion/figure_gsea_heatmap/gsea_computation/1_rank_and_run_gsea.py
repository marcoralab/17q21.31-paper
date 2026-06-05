import os
import pandas as pd
import gseapy as gp
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# List of all input files
input_files = [
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1A_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1A_vs_AA_h1h1_mfg/all_genes_result_cov_rin_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1A_vs_AA_h1h1_caudate/all_genes_result_cov_rin_sex_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1A_vs_AA_h1h1_caudate/all_genes_result_cov_rin_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_AA_h1h1_mfg/all_genes_result_cov_rin_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_AA_h1h1_caudate/all_genes_result_cov_rin_hap_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_AA_h1h1_caudate/all_genes_result_cov_rin_sex_hap_with_name_and_type.csv",
    "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_mfg/all_genes_result_cov_rin_ancestry_with_name_and_type.csv",
    "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_caudate_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type.csv",
    "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_caudate/all_genes_result_cov_rin_ancestry_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_sex_ancestry_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_ancestry_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_caudate/all_genes_result_cov_rin_sex_ancestry_with_name_and_type.csv",
 #   "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_caudate/all_genes_result_cov_rin_ancestry_with_name_and_type.csv",
    "path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/EUR_h1h2_vs_EUR_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type.csv"
]

def process_file(file_path):
    try:
        # Extract the parent folder and file name for organization
        parent_folder = Path(file_path).parts[-2]
        file_name = Path(file_path).stem
        
        # Create output directory
        output_dir = f"GSEA_Analysis_{parent_folder}/{file_name}"
        os.makedirs(output_dir, exist_ok=True)
        
        # Create subdirectories for different GSEA results
        go_bp_dir = os.path.join(output_dir, 'GSEA_results_GO_BP')
        kegg_dir = os.path.join(output_dir, 'GSEA_results_KEGG')
        reactome_dir = os.path.join(output_dir, 'GSEA_results_reactome')
        
        # Load and process data
        df = pd.read_csv(file_path)
        
        # Log initial gene count
        initial_count = len(df)
        logger.info(f"Initial gene count: {initial_count}")

        # Track each filtering step separately
        na_pvalue_count = df['pvalue'].isna().sum()
        na_stat_count = df['stat'].isna().sum()
        na_log2fc_count = df['log2FoldChange'].isna().sum()

        # Remove genes with NA values in key columns
        df = df.dropna(subset=['pvalue', 'stat', 'log2FoldChange'])
        after_na_count = len(df)
        na_removed_count = initial_count - after_na_count 

        logger.info(f"Removed {na_pvalue_count} genes with NA p-values")
        logger.info(f"Removed {na_stat_count} genes with NA statistics")
        logger.info(f"Removed {na_log2fc_count} genes with NA log2FoldChange")
        logger.info(f"Total genes after NA filtering: {after_na_count}")

        # Convert Gene Symbols to uppercase and remove any whitespace
        df['Gene_Symbol'] = df['Gene_Symbol'].str.upper().str.strip()
        
        # Remove any genes with missing gene symbols
        before_symbol_filter = len(df)
        df = df[df['Gene_Symbol'].notna()]
        df = df[df['Gene_Symbol'] != '']
        symbol_removed_count = before_symbol_filter - len(df)
        logger.info(f"Removed {symbol_removed_count} genes with missing symbols")
        
        # Use DESeq2 statistic for ranking
        df['rank_metric'] = df['stat']
        df['rank'] = df['rank_metric'].rank(ascending=False)
        ranked_genes = df.sort_values('rank', ascending=True)
        
        # Create ranked gene list
        ranked_gene_list = ranked_genes[['Gene_Symbol', 'rank_metric']].set_index('Gene_Symbol').squeeze()
        ranked_list_path = os.path.join(output_dir, "ranked_gene_list_for_gsea.rnk")
        ranked_gene_list.to_csv(ranked_list_path, sep="\t", header=False, index=True)
        
        # Skip GSEA for inversion_plus_1MB files
        if 'inversion_plus_1MB' in file_path:
            logger.info(f"Skipping GSEA for inversion file: {file_name}")
            
            # Still save the ranked genes and statistics
            ranked_genes[['Gene_Symbol', 'rank_metric', 'pvalue', 'log2FoldChange']].to_csv(
                os.path.join(output_dir, "ranked_genes.csv"), index=False
            )
            return True
        
        # Common GSEA parameters
        gsea_kwargs = {
            'rnk': ranked_list_path,
            'permutation_num': 1000,  # Increased for better statistics
            'min_size': 15,   # Standard minimum from MSigDB
            'max_size': 500,  # Standard maximum from MSigDB
            'processes': 4,   # Use multiple cores
            'seed': 42,      # For reproducibility
            'no_plot': True  # Skip plotting to save time
        }
        
        # Try GO BP analysis
        try:
            gsea_results_go_bp = gp.prerank(
                gene_sets='path/to/msigdb/c5.go.bp.v2023.2.Hs.symbols.gmt',
                outdir=go_bp_dir,
                **gsea_kwargs
            )
            logger.info(f"Completed GO BP analysis for {file_name}")
        except Exception as e:
            logger.error(f"Error in GO BP analysis for {file_name}: {str(e)}")
        
        # Try KEGG analysis
        try:
            gsea_results_kegg = gp.prerank(
                gene_sets='path/to/msigdb/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt',
                outdir=kegg_dir,
                **gsea_kwargs
            )
            logger.info(f"Completed KEGG analysis for {file_name}")
        except Exception as e:
            logger.error(f"Error in KEGG analysis for {file_name}: {str(e)}")
        
        # Try Reactome analysis
        try:
            gsea_results_reactome = gp.prerank(
                gene_sets='path/to/msigdb/c2.cp.reactome.v2023.2.Hs.symbols.gmt',
                outdir=reactome_dir,
                **gsea_kwargs
            )
            logger.info(f"Completed Reactome analysis for {file_name}")
        except Exception as e:
            logger.error(f"Error in Reactome analysis for {file_name}: {str(e)}")
        
        # Save all ranked genes with full information
        ranked_genes[['Gene_Symbol', 'rank_metric', 'pvalue', 'padj', 'log2FoldChange']].to_csv(
            os.path.join(output_dir, "ranked_genes.csv"), index=False
        )
        
        # Save top ranked genes
        top_genes = ranked_genes[['Gene_Symbol', 'rank_metric', 'pvalue', 'log2FoldChange']].head(50)
        top_genes.to_csv(os.path.join(output_dir, "top_50_ranked_genes.csv"), index=False)
        
        # Create summary report
        with open(os.path.join(output_dir, "analysis_summary.txt"), 'w') as f:
            f.write(f"GSEA Analysis Summary for {parent_folder}/{file_name}\n")
            f.write("=" * 50 + "\n\n")
            f.write("Analysis completed on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
            f.write(f"Input file: {file_path}\n")
            f.write(f"Initial number of genes: {initial_count}\n")
            f.write(f"Genes removed due to missing values: {na_removed_count}\n")
            f.write(f"Genes removed due to missing symbols: {symbol_removed_count}\n")
            f.write(f"Final number of genes analyzed: {len(ranked_genes)}\n")
            f.write(f"Statistic range: {ranked_genes['rank_metric'].min():.3f} to {ranked_genes['rank_metric'].max():.3f}\n")
            f.write("\nGene set sizes:\n")
            f.write(f"Minimum: {gsea_kwargs['min_size']}\n")
            f.write(f"Maximum: {gsea_kwargs['max_size']}\n")
        
        logger.info(f"Successfully processed {parent_folder}/{file_name}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing {file_path}: {str(e)}")
        return False

def main():
    # Create main output directory
    os.makedirs("GSEA_Analyses", exist_ok=True)
    
    # Process all files in parallel
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(process_file, input_files))
    
    # Create overall summary
    successful = sum(results)
    failed = len(input_files) - successful
    
    with open("GSEA_Analyses/overall_summary.txt", 'w') as f:
        f.write("GSEA Analyses Overall Summary\n")
        f.write("=" * 30 + "\n\n")
        f.write(f"Total files processed: {len(input_files)}\n")
        f.write(f"Successful: {successful}\n")
        f.write(f"Failed: {failed}\n")
        f.write("\nProcessed files:\n")
        for file_path in input_files:
            parent = Path(file_path).parts[-2]
            name = Path(file_path).stem
            f.write(f"- {parent}/{name}\n")

if __name__ == "__main__":
    main()