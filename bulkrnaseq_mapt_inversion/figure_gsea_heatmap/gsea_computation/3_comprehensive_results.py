import pandas as pd
import os
from typing import Dict, Set, List

def read_gsea_results(filepath):
    """Read GSEA results from CSV file and return dataframe with all columns."""
    df = pd.read_csv(filepath)
    return df

def find_common_pathways(results_dict, conditions, direction='up'):
    """Find pathways that meet specified direction criteria across conditions."""
    first_cond = conditions[0]
    if direction == 'up':
        pathways = set(results_dict[first_cond][results_dict[first_cond]['NES'] > 0]['Term'])
    else:
        pathways = set(results_dict[first_cond][results_dict[first_cond]['NES'] < 0]['Term'])
    
    for cond in conditions[1:]:
        if direction == 'up':
            curr_pathways = set(results_dict[cond][results_dict[cond]['NES'] > 0]['Term'])
        else:
            curr_pathways = set(results_dict[cond][results_dict[cond]['NES'] < 0]['Term'])
        pathways = pathways.intersection(curr_pathways)
    
    return pathways

def find_opposite_pathways(results_dict, up_conditions, down_conditions):
    """Find pathways that are up in some conditions and down in others."""
    up_pathways = set()
    for cond in up_conditions:
        curr_pathways = set(results_dict[cond][results_dict[cond]['NES'] > 0]['Term'])
        if not up_pathways:
            up_pathways = curr_pathways
        else:
            up_pathways = up_pathways.intersection(curr_pathways)
    
    down_pathways = set()
    for cond in down_conditions:
        curr_pathways = set(results_dict[cond][results_dict[cond]['NES'] < 0]['Term'])
        if not down_pathways:
            down_pathways = curr_pathways
        else:
            down_pathways = down_pathways.intersection(curr_pathways)
    
    return up_pathways.intersection(down_pathways)

def find_mixed_effects(pathway_sets: Dict[str, Set]):
    """Find pathways that appear in multiple criteria groups."""
    all_pathways = set().union(*pathway_sets.values())
    mixed_effects = {}
    for pathway in all_pathways:
        patterns = []
        for set_name, pathway_set in pathway_sets.items():
            if pathway in pathway_set:
                patterns.append(set_name)
        if len(patterns) > 1:
            mixed_effects[pathway] = patterns
    return mixed_effects

def get_pathway_genes(results_dict: Dict, pathway: str, condition: str) -> str:
    """Extract genes for a specific pathway and condition."""
    df = results_dict[condition]
    pathway_row = df[df['Term'] == pathway]
    if not pathway_row.empty and 'Lead_genes' in pathway_row.columns:
        return pathway_row['Lead_genes'].iloc[0]
    return ''

def create_comprehensive_table(results_dict: Dict, pathway_sets: Dict, mixed_effects: Dict):
    """Create a comprehensive table with pathway information."""
    rows = []
    conditions = list(results_dict.keys())
    
    # Process each pathway set
    for effect_type, pathways in pathway_sets.items():
        for pathway in pathways:
            row = {
                'Pathway': pathway,
                'Effect_Type': effect_type,
                'Mixed_Effects': 'No',
                'Mixed_Effect_Pattern': '',
                'NES_Pattern': ''
            }
            
            # Add genes and NES for each condition
            for condition in conditions:
                df = results_dict[condition]
                pathway_row = df[df['Term'] == pathway]
                
                if not pathway_row.empty:
                    row[f'{condition}_genes'] = pathway_row['Lead_genes'].iloc[0] if 'Lead_genes' in pathway_row.columns else ''
                    row[f'{condition}_NES'] = pathway_row['NES'].iloc[0]
                    row[f'{condition}_FDR'] = pathway_row['FDR q-val'].iloc[0]
                else:
                    row[f'{condition}_genes'] = ''
                    row[f'{condition}_NES'] = ''
                    row[f'{condition}_FDR'] = ''
            
            # Add mixed effects information if applicable
            if pathway in mixed_effects:
                row['Mixed_Effects'] = 'Yes'
                row['Mixed_Effect_Pattern'] = '; '.join(mixed_effects[pathway])
            
            # Create NES pattern description
            nes_patterns = []
            for condition in conditions:
                if row[f'{condition}_NES']:
                    direction = 'UP' if float(row[f'{condition}_NES']) > 0 else 'DOWN'
                    nes_patterns.append(f"{condition}: {direction}")
            row['NES_Pattern'] = '; '.join(nes_patterns)
            
            rows.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(rows)
    
    # Reorder columns
    column_order = ['Pathway', 'Effect_Type', 'Mixed_Effects', 'Mixed_Effect_Pattern', 'NES_Pattern']
    for condition in conditions:
        column_order.extend([f'{condition}_NES', f'{condition}_FDR', f'{condition}_genes'])
    
    df = df[column_order]
    
    return df

def main():
    # Define file paths
   # result_files = {
   #     "AA_H1H1A_vs_H1H1": "neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_A_aa_vs_h1h1_aa/all_genes_results_h1h1_A_aa_vs_h1h1_aa_cov_haplotype_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H2_vs_H1H1": "neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h1_aa/all_genes_results_h1h2_vs_h1h1_aa_cov_haplotype_sex_rin_with_gene_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H1_vs_EUR_H1H1": "neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_aa_vs_h1h1_eur/all_genes_results_h1h1_aa_vs_eur_cov_ancestry_sex_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H2_vs_EUR_H1H2": "neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h2_eur/all_genes_results_cov_haplotype_sex_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "EUR_H1H2_vs_H1H1": "neurons/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_vs_h1h1_eur/all_genes_results_cov_haplotype_sex_rin_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv"
   # }

   # BRAIN:
   # result_files = {
   #     "AA_H1H1A_vs_H1H1": "brain/no_na_padj/GSEA_Analysis_AA_h1h1A_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H2_vs_H1H1": "brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_AA_h1h1_mfg/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H1_vs_EUR_H1H1": "brain/no_na_padj/GSEA_Analysis_AA_h1h1_vs_eur_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H2_vs_EUR_H1H2": "brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_eur_h1h2_mfg/all_genes_result_cov_rin_sex_ancestry_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "EUR_H1H2_vs_H1H1": "brain/no_na_padj/GSEA_Analysis_EUR_h1h2_vs_EUR_h1h1_mfg_no331/all_genes_result_cov_rin_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv"
   # }
        
   #ASTRO
   # result_files = {
   #     "AA_H1H1A_vs_H1H1": "astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_and_h1h1a_aa/all_genes_results_H1H1_A_vs_H1H1_AA_cov_haplotype_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H2_vs_H1H1": "astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_aa_vs_h1h1_aa/all_genes_results_H1H2_vs_H1H1_AA_cov_haplotype_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H1_vs_EUR_H1H1": "astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h1_eur_and_aa/all_genes_results_AAvsEURH1H1_cov_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "AA_H1H2_vs_EUR_H1H2": "astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_and_aa/all_genes_results_h1h2_aa_vs_eur_cov_ancestry_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
   #     "EUR_H1H2_vs_H1H1": "astro/only_covariates_excel_number_1_option/no_na_padj/GSEA_Analysis_h1h2_eur_vs_h1h1_eur/all_genes_results_H1H2_vs_H1H1_EUR_cov_haplotype_sex_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv"
 #Brain caudate
    result_files = {
        "AA_H1H1A_vs_H1H1": "brain/no_na_padj/GSEA_Analysis_AA_h1h1A_vs_AA_h1h1_caudate/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
        "AA_H1H2_vs_H1H1": "brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_AA_h1h1_caudate/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
        "EUR_H1H1_vs_AA_H1H1": "brain/no_na_padj/GSEA_Analysis_AA_h1h1_vs_eur_h1h1_caudate/all_genes_result_cov_rin_ancestry_sex_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
        "EUR_H1H2_vs_AA_H1H2": "brain/no_na_padj/GSEA_Analysis_AA_h1h2_vs_eur_h1h2_caudate/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv",
        "EUR_H1H2_vs_H1H1": "brain/no_na_padj/GSEA_Analysis_EUR_h1h2_vs_EUR_h1h1_caudate/all_genes_result_cov_rin_sex_hap_with_name_and_type/GSEA_results_KEGG_no_na_padj/gseapy.gene_set.prerank.report.csv"
    
    } 
        
    # Read all results
    results_dict = {}
    for name, filepath in result_files.items():
        if os.path.exists(filepath):
            results_dict[name] = read_gsea_results(filepath)
        else:
            print(f"Warning: File not found - {filepath}")
    
    pathway_sets = {}
    
    # Find pathways for different conditions
    up_conditions_three = ["AA_H1H1A_vs_H1H1", "AA_H1H2_vs_H1H1", "EUR_H1H2_vs_H1H1"]
    pathway_sets['up_807SNP_effect'] = find_common_pathways(results_dict, up_conditions_three, 'up')
    pathway_sets['down_807SNP_effect'] = find_common_pathways(results_dict, up_conditions_three, 'down')
     
    up_conditions_h1h2 = ["AA_H1H2_vs_H1H1"]
    down_conditions_aa_eur = ["AA_H1H1A_vs_H1H1", "EUR_H1H2_vs_H1H1", "EUR_H1H2_vs_AA_H1H2"]
    pathway_sets['up_H2_in_AFR_background'] = find_opposite_pathways(results_dict, up_conditions_h1h2, down_conditions_aa_eur)
    pathway_sets['down_H2_in_AFR_background'] = find_opposite_pathways(results_dict, down_conditions_aa_eur, up_conditions_h1h2)
    
    up_conditions_new = ["AA_H1H2_vs_H1H1", "EUR_H1H2_vs_H1H1"]
    down_conditions_new = ["AA_H1H1A_vs_H1H1"]
    pathway_sets['up_haplotype_effect'] = find_opposite_pathways(results_dict, up_conditions_new, down_conditions_new)
    pathway_sets['down_haplotype_effect'] = find_opposite_pathways(results_dict, down_conditions_new, up_conditions_new)
    
    # Find mixed effects
    mixed_effects = find_mixed_effects(pathway_sets)
    
    # Create comprehensive table
    results_table = create_comprehensive_table(results_dict, pathway_sets, mixed_effects)
    
    # Save results
    results_table.to_csv('gsea_comprehensive_results_Brain_caudate_no_na_padj_flipped_last2.csv', index=False)
    print("Results saved to 'gsea_comprehensive_results_cau.csv'")

if __name__ == "__main__":
    main()