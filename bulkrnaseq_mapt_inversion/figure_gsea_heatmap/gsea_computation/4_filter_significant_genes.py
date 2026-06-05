import pandas as pd
import os
from itertools import product
import re
from collections import defaultdict

def extract_gene_names(genes_str):
    """Extract just the gene names from the genes string."""
    if pd.isna(genes_str) or genes_str == "Not found in this comparison" or genes_str == "No genes with p-value < 0.05":
        return set()
    
    genes = genes_str.split("; ")
    gene_names = set()
    for gene in genes:
        match = re.match(r'(\w+)\s+\(', gene)
        if match:
            gene_names.add(match.group(1))
    return gene_names

def find_shared_genes(df, pathway, min_datasets=2):
    """Find genes that appear in at least min_datasets for a given pathway."""
    # Get all rows for this pathway
    pathway_rows = df[df['Pathway'] == pathway]
    
    # Create a dict to count gene occurrences
    gene_counts = defaultdict(int)
    
    # Count occurrences of each gene
    for _, row in pathway_rows.iterrows():
        genes = extract_gene_names(row['Genes'])
        for gene in genes:
            gene_counts[gene] += 1
    
    # Get genes that appear in at least min_datasets files
    shared_genes = {gene for gene, count in gene_counts.items() if count >= min_datasets}
    
    return shared_genes if shared_genes else set()

def add_shared_genes_columns(df):
    """Add columns for genes shared across different numbers of datasets."""
    pathways = df['Pathway'].unique()
    
    # Initialize new columns
    df['shared_genes_across_5'] = ''
    df['shared_genes_across_4'] = ''
    df['shared_genes_across_3'] = ''
    df['shared_genes_across_2'] = ''
    
    for pathway in pathways:
        shared_5 = find_shared_genes(df, pathway, 5)
        shared_4 = find_shared_genes(df, pathway, 4)
        shared_3 = find_shared_genes(df, pathway, 3)
        shared_2 = find_shared_genes(df, pathway, 2)
        
        # Remove overlapping genes
        shared_4 = shared_4 - shared_5
        shared_3 = shared_3 - shared_4 - shared_5
        shared_2 = shared_2 - shared_3 - shared_4 - shared_5
        
        # Update all rows for this pathway
        df.loc[df['Pathway'] == pathway, 'shared_genes_across_5'] = '; '.join(sorted(shared_5)) if shared_5 else ''
        df.loc[df['Pathway'] == pathway, 'shared_genes_across_4'] = '; '.join(sorted(shared_4)) if shared_4 else ''
        df.loc[df['Pathway'] == pathway, 'shared_genes_across_3'] = '; '.join(sorted(shared_3)) if shared_3 else ''
        df.loc[df['Pathway'] == pathway, 'shared_genes_across_2'] = '; '.join(sorted(shared_2)) if shared_2 else ''
    
    return df

def filter_significant_genes(genes_str, pvalue_threshold=0.05):
    """Filter genes to only include those with p-value < threshold."""
    if pd.isna(genes_str) or genes_str == "Not found in this comparison":
        return genes_str
        
    # Split genes string into individual genes
    genes = genes_str.split("; ")
    significant_genes = []
    
    for gene in genes:
        # Extract p-value using regex
        match = re.search(r'p-value: ([\d.e-]+)', gene)
        if match:
            pvalue = float(match.group(1))
            if pvalue < pvalue_threshold:
                significant_genes.append(gene)
    
    if significant_genes:
        return "; ".join(significant_genes)
    else:
        return "No genes with p-value < 0.05"

def get_pathways_from_comprehensive_results(comprehensive_file):
    """Extract unique pathways from a comprehensive results file."""
    try:
        df = pd.read_csv(comprehensive_file)
        return df['Pathway'].unique().tolist()
    except Exception as e:
        print(f"Error reading {comprehensive_file}: {str(e)}")
        return []

def process_tissue_files(tissue_files, pathways_of_interest, base_path):
    """Process all files for a tissue type and ensure all pathway-file combinations exist."""
    all_results = []
    
    # Use just filenames for the template
    file_basenames = [os.path.basename(f) for f in tissue_files]
    combinations = list(product(pathways_of_interest, file_basenames))
    template_df = pd.DataFrame(combinations, columns=['Pathway', 'Source_File'])
    
    # Process each file
    for file in tissue_files:
        if not os.path.exists(file):
            print(f"Warning: File {file} not found, skipping...")
            continue
            
        try:
            df = pd.read_csv(file)
            file_basename = os.path.basename(file)
            
            # Filter for pathways of interest
            filtered_df = df[df['Pathway'].isin(pathways_of_interest)].copy()
            
            if not filtered_df.empty:
                # Filter genes by p-value
                filtered_df['Genes'] = filtered_df['Genes'].apply(filter_significant_genes)
                
                filtered_df['Source_File'] = file_basename
                all_results.append(filtered_df)
            
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
    
    # Combine all results
    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)
    else:
        results_df = pd.DataFrame(columns=['Pathway', 'NES', 'FDR q-value', 'Genes', 'Source_File'])
    
    # Merge with template to ensure all combinations exist
    complete_df = template_df.merge(results_df, on=['Pathway', 'Source_File'], how='left')
    
    # Fill NaN values appropriately
    complete_df['NES'] = complete_df['NES'].fillna(0)
    complete_df['FDR q-value'] = complete_df['FDR q-value'].fillna(1)
    complete_df['Genes'] = complete_df['Genes'].fillna("Not found in this comparison")
    
    # Add shared genes columns
    complete_df = add_shared_genes_columns(complete_df)
    
    return complete_df

def main():
    base_path = "path/to/BulkRNASeq/Brain_data/gsea_caudate_mfg_split/pathway_analysis_results"
    
    # Comprehensive results files for each tissue
    comprehensive_files = {
        'neurons': 'gsea_comprehensive_results_neurons_no_na_padj_flipped_last2.csv',
      #  'mfg': 'gsea_comprehensive_results_Brain_MFG_no_na_padj_flipped_last2.csv',
        'mfg': 'gsea_comprehensive_results_mfg_added_2_samples_no_na_padj_flipped_last2.csv',
        'astro': 'gsea_comprehensive_results_astro_no_na_padj_flipped_last2.csv',
        'caudate': 'gsea_comprehensive_results_Brain_caudate_no_na_padj_flipped_last2.csv'
    }

    # File patterns for each tissue type
    file_patterns = {
    
          'mfg': [
            os.path.join(base_path, 'eur_h1h2_vs_eur_h1h1_brain_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),  # unchanged
            os.path.join(base_path, 'aa_h1h2_vs_aa_h1h1_brain_2_added_samples_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h1a_vs_aa_h1h1_brain_2_added_samples_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa_2_added_samples_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa_2_added_samples_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv')
        ],

     #   'mfg': [
     #       os.path.join(base_path, 'eur_h1h2_vs_eur_h1h1_brain_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
     #       os.path.join(base_path, 'aa_h1h2_vs_aa_h1h1_brain_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
     #       os.path.join(base_path, 'aa_h1h1a_vs_aa_h1h1_brain_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
     #       os.path.join(base_path, 'h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
     #       os.path.join(base_path, 'h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv')
     #   ],
        'neurons': [
            os.path.join(base_path, 'eur_h1h2_vs_eur_h1h1_neuron_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h2_vs_aa_h1h1_neuron_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h1a_vs_aa_h1h1_neuron_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h1_aa_vs_h1h1_eur_neurons_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h2_aa_vs_h1h2_eur_neurons_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv')
        ],
        'astro': [
            os.path.join(base_path, 'eur_h1h2_vs_eur_h1h1_astro_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h2_vs_aa_h1h1_astro_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h1a_vs_aa_h1h1_astro_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h1_aa_vs_h1h1_eur_astro_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h2_aa_vs_h1h2_eur_astro_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv')
        ],
        'caudate': [
            os.path.join(base_path, 'h1h2_aa_vs_h1h2_eur_brain_caudate_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'h1h1_aa_vs_h1h1_eur_brain_caudate_flipped_eur_vs_aa_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'eur_h1h2_vs_eur_h1h1_brain_caudate_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h2_vs_aa_h1h1_brain_caudate_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv'),
            os.path.join(base_path, 'aa_h1h1a_vs_aa_h1h1_brain_caudate_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv')
        ]
    }

    # Process each tissue type using pathways from comprehensive results
    for tissue, comp_file in comprehensive_files.items():
        print(f"\nProcessing {tissue} files...")
        
        # Get pathways from comprehensive results file
        pathways = get_pathways_from_comprehensive_results(comp_file)
        if not pathways:
            print(f"No pathways found in {comp_file}, skipping tissue...")
            continue
            
        print(f"Found {len(pathways)} pathways in comprehensive results")
        
        # Process files with these pathways
        results_df = process_tissue_files(file_patterns[tissue], pathways, base_path)
        
        if not results_df.empty:
            output_filename = f'{tissue}_filtered_pathways_sig_genes_from_comprehensive_gsea_added_2_mfg_samples.csv'
            results_df.to_csv(output_filename, index=False)
            print(f"Created {output_filename} with {len(results_df)} rows")
            # Print summary of pathways found
            pathway_counts = results_df.groupby('Pathway').size()
            print("\nPathway counts:")
            print(pathway_counts)
        else:
            print(f"No matching pathways found for {tissue}")

if __name__ == "__main__":
    main()