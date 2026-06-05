import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def format_small_float(value):
    """
    Format float values to show actual values in scientific notation,
    no matter how small they are
    """
    try:
        value = float(value)
        if value == 0:  # Only if it's exactly zero
            return "0.00e+00"
        else:
            return '{:.2e}'.format(value)
    except (ValueError, TypeError):
        return "N/A"

def generate_csv(label, gsea_path, gene_path, output_csv):
    # Load GSEA results and gene-level data
    gsea_df = pd.read_csv(gsea_path)
    gene_data = pd.read_csv(gene_path)

    # Initialize a list to store the output rows
    output_rows = []

    # Iterate over GSEA result rows
    for _, gsea_row in gsea_df.iterrows():
        pathway_name = gsea_row['Term']
        nes = gsea_row['NES']
        
        # Format FDR q-value showing actual values no matter how small
        fdr = format_small_float(gsea_row['FDR q-val'])
        
        genes = gsea_row['Lead_genes'].split(';')

        # For each gene, get log2FoldChange, p-value, and padj from the gene data
        gene_info = []
        for gene in genes:
            gene = gene.strip()  # Remove any whitespace
            log2fc_value = gene_data.loc[gene_data['Gene_Symbol'] == gene, 'log2FoldChange'].values
            pvalue_value = gene_data.loc[gene_data['Gene_Symbol'] == gene, 'pvalue'].values
            padj_value = gene_data.loc[gene_data['Gene_Symbol'] == gene, 'padj'].values

            if log2fc_value.size > 0 and pvalue_value.size > 0 and padj_value.size > 0:
                gene_info.append(f"{gene} (log2FC: {log2fc_value[0]:.4f}, p-value: {format_small_float(pvalue_value[0])}, adj p-value: {format_small_float(padj_value[0])})")
            elif log2fc_value.size > 0 and pvalue_value.size > 0:
                gene_info.append(f"{gene} (log2FC: {log2fc_value[0]:.4f}, p-value: {format_small_float(pvalue_value[0])}, adj p-value: N/A)")
            elif log2fc_value.size > 0:
                gene_info.append(f"{gene} (log2FC: {log2fc_value[0]:.4f}, p-value: N/A, adj p-value: N/A)")
            else:
                gene_info.append(f"{gene} (log2FC: N/A, p-value: N/A, adj p-value: N/A)")

        # Join all gene information into one string
        gene_info_str = "; ".join(gene_info)

        # Append the information as a new row to the output
        output_rows.append({
            'Pathway': pathway_name,
            'NES': '{:.4f}'.format(nes),
            'FDR q-value': fdr,
            'Genes': gene_info_str
        })

    # Convert output rows to a DataFrame and save as CSV
    output_df = pd.DataFrame(output_rows)
    # Sort by absolute NES value
    output_df = output_df.sort_values(by='NES', key=lambda x: abs(pd.to_numeric(x)), ascending=False)
    output_df.to_csv(output_csv, index=False)
    print(f"CSV saved to {output_csv}")

# List of GSEA result files and corresponding gene-level data

gsea_files = {
      'aa_h1h1a_vs_aa_h1h1_brain_2_added_samples': (
          'brain/no_na_padj/GSEA_Analysis_AA_h1h1A_vs_AA_h1h1_mfg/2_added_samples_Ensembl111/all_genes_result_cov_rin_hap_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv',
          'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1A_vs_AA_h1h1_mfg/2_added_samples_Ensembl111/all_genes_result_cov_rin_hap_sex_with_name_and_type.csv'
      ),
      'aa_h1h2_vs_aa_h1h1_brain_2_added_samples': (
          'brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_AA_h1h1_mfg/2_added_samples_Ensembl111/all_genes_result_cov_rin_hap_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv',
          'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_AA_h1h1_mfg/2_added_samples_Ensembl111/all_genes_result_cov_rin_hap_sex_with_name_and_type.csv'
      ),
      'h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa_2_added_samples': (
          'brain/no_na_padj/GSEA_Analysis_AA_h1h1_vs_eur_h1h1_mfg_no331/2_added_samples_Ensembl111/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv',
          'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_mfg_no331/2_added_samples_Ensembl111/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type.csv'
      ),
      'h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa_2_added_samples': (
          'brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_eur_h1h2_mfg/2_added_samples_Ensembl111/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv',
          'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_mfg/2_added_samples_Ensembl111/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type.csv'
      ),



#gsea_files = {
#    #CAUDATE:
    'h1h1_aa_vs_h1h1_eur_brain_caudate_flipped_eur_vs_aa': ('brain/no_na_padj/GSEA_Analysis_AA_h1h1_vs_eur_h1h1_caudate/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_caudate/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type.csv'),
    'h1h2_aa_vs_h1h2_eur_brain_caudate_flipped_eur_vs_aa': ('brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_eur_h1h2_caudate/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_caudate/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type.csv'),
    
   # 'h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa': ('brain/no_na_padj/GSEA_Analysis_AA_h1h1_vs_eur_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv','path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type.csv'),
   # 'h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa': ('brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type.csv'),
    
    'h1h1_aa_vs_h1h1_eur_neurons_flipped_eur_vs_aa': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_aa_vs_h1h1_eur/all_genes_results_h1h1_aa_vs_eur_cov_ancestry_sex_rin_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h1_aa_vs_h1h1_eur/all_genes_results_h1h1_aa_vs_eur_cov_ancestry_sex_rin_flipped_eur_vs_aa_with_name_and_type.csv'),
    'h1h2_aa_vs_h1h2_eur_neurons_flipped_eur_vs_aa': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h2_eur/all_genes_results_cov_haplotype_sex_rin_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h2_aa_vs_h1h2_eur/all_genes_results_cov_haplotype_sex_rin_flipped_eur_vs_aa_with_name_and_type.csv'),


    'h1h1_aa_vs_h1h1_eur_astro_flipped_eur_vs_aa': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_eur_and_aa/all_genes_results_AAvsEURH1H1_cov_ancestry_sex_flipped_eur_vs_aa_with_name_and_type_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h1_eur_and_aa/all_genes_results_AAvsEURH1H1_cov_ancestry_sex_flipped_eur_vs_aa_with_name_and_type_with_name_and_type.csv'),
    'h1h2_aa_vs_h1h2_eur_astro_flipped_eur_vs_aa': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_and_aa/all_genes_results_h1h2_aa_vs_eur_cov_ancestry_sex_flipped_eur_vs_aa_with_name_and_type_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h2_eur_and_aa/all_genes_results_h1h2_aa_vs_eur_cov_ancestry_sex_flipped_eur_vs_aa_with_name_and_type_with_name_and_type.csv'),

    
   # 'aa_h1h1a_vs_aa_h1h1_brain': ('brain/no_na_padj/GSEA_Analysis_AA_h1h1A_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1A_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type.csv'),
   # 'aa_h1h2_vs_aa_h1h1_brain': ('brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type.csv'),
    'eur_h1h2_vs_eur_h1h1_brain': ('brain/no_na_padj/GSEA_Analysis_EUR_h1h2_vs_EUR_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/EUR_h1h2_vs_EUR_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type.csv'),
    'aa_h1h1a_vs_aa_h1h1_astro': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_and_h1h1a_aa/all_genes_results_H1H1_A_vs_H1H1_AA_cov_haplotype_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h1_and_h1h1a_aa/all_genes_results_H1H1_A_vs_H1H1_AA_cov_haplotype_with_name_and_type.csv'),
    'aa_h1h2_vs_aa_h1h1_astro': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h1_aa/all_genes_results_H1H2_vs_H1H1_AA_cov_haplotype_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h2_aa_vs_h1h1_aa/all_genes_results_H1H2_vs_H1H1_AA_cov_haplotype_sex_with_name_and_type.csv'),
    'eur_h1h2_vs_eur_h1h1_astro': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_vs_h1h1_eur/all_genes_results_H1H2_vs_H1H1_EUR_cov_haplotype_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h2_eur_vs_h1h1_eur/all_genes_results_H1H2_vs_H1H1_EUR_cov_haplotype_sex_with_name_and_type.csv'),
  
  #  'h1h1_aa_vs_h1h1_eur_brain': ('brain/no_na_padj/GSEA_Analysis_AA_h1h1_vs_eur_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h1_vs_eur_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type.csv'),
   # 'h1h2_aa_vs_h1h2_eur_brain': ('brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_sex_ancestry_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_sex_ancestry_with_name_and_type.csv'),

  #  'h1h1_aa_vs_h1h1_eur_astro': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_eur_and_aa/all_genes_results_AAvsEURH1H1_cov_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h1_eur_and_aa/all_genes_results_AAvsEURH1H1_cov_ancestry_sex_with_name_and_type.csv'),
  #  'h1h2_aa_vs_h1h2_eur_astro': ('astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_and_aa/all_genes_results_h1h2_aa_vs_eur_cov_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/astro_re/Routputs_haplotype_specific/rsem_redo/h1h2_eur_and_aa/all_genes_results_h1h2_aa_vs_eur_cov_ancestry_sex_with_name_and_type.csv'),

  #  'h1h1_aa_vs_h1h1_eur_neurons': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_aa_vs_h1h1_eur/all_genes_results_h1h1_aa_vs_eur_cov_ancestry_sex_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h1_aa_vs_h1h1_eur/all_genes_results_h1h1_aa_vs_eur_cov_ancestry_sex_rin_with_name_and_type.csv'),
  #  'h1h2_aa_vs_h1h2_eur_neurons': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h2_eur/all_genes_results_cov_haplotype_sex_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h2_aa_vs_h1h2_eur/all_genes_results_cov_haplotype_sex_rin_with_name_and_type.csv')
    'aa_h1h2_vs_aa_h1h1_neuron': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h1_aa/all_genes_results_h1h2_vs_h1h1_aa_cov_haplotype_sex_rin_with_gene_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h2_aa_vs_h1h1_aa/all_genes_results_h1h2_vs_h1h1_aa_cov_haplotype_sex_rin_with_gene_name_and_type.csv'),
    'eur_h1h2_vs_eur_h1h1_neuron': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_vs_h1h1_eur/all_genes_results_cov_haplotype_sex_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h2_eur_vs_h1h1_eur/all_genes_results_cov_haplotype_sex_rin_with_name_and_type.csv'),
      'aa_h1h1a_vs_aa_h1h1_neuron': ('neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_A_aa_vs_h1h1_aa/all_genes_results_h1h1_A_aa_vs_h1h1_aa_cov_haplotype_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv', 'path/to/BulkRNASeq/Ancestry_analysis/neurons_re/Routputs_haplotype_specific/rsem_redo/h1h1_A_aa_vs_h1h1_aa/all_genes_results_h1h1_A_aa_vs_h1h1_aa_cov_haplotype_rin_with_name_and_type.csv'),

}

# Loop through each GSEA result and generate a separate CSV file for each dataset
for label, (gsea_path, gene_path) in gsea_files.items():
    output_csv = f"pathway_analysis_results/rerun_again_confirmation/{label}_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv"
    generate_csv(label, gsea_path, gene_path, output_csv)