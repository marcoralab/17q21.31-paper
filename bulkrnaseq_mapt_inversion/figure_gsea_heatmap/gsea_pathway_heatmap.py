import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# ===============================================================
# INPUTS / OUTPUTS
# This script reads its GSEA result CSVs by filename from the current working
# directory, and writes the heatmap PNGs there too. Run it from the directory
# that contains the GSEA result files. Required inputs (per tissue):
#   - gsea_comprehensive_results_<tissue>_no_na_padj_flipped_last2.csv
#       (used for pathway ordering; for MFG the "_added_2_samples" variant)
#   - <tissue>_filtered_pathways_sig_genes_from_comprehensive_gsea[...].csv
#       (the per-tissue significant-pathway gene table that is plotted)
# The paper figure is `mfg_heatmap_2_added_samples.png` (tissue_name='MFG').
# ===============================================================

def clean_pathway_name(name):
    if name.startswith('KEGG_'):
        name = name[5:]
    return name.replace('_', ' ')

def get_comparison_type(filename):
    # Standardize the comparison names
    if 'h1h1a_vs_h1h1' in filename or 'aa_h1h1a_vs_aa_h1h1' in filename:
        return 'AA H1H1A vs AA H1H1'
    elif 'eur_h1h2_vs_eur_h1h1' in filename:
        return 'EUR H1H2 vs EUR H1H1'
    elif 'aa_h1h2_vs_aa_h1h1' in filename:
        return 'AA H1H2 vs AA H1H1'
    elif 'h1h1_aa_vs_h1h1_eur' in filename or 'h1h1_aa_vs_h1h1_eur_brain' in filename:
        return 'H1H1 AA vs H1H1 EUR'
    elif 'h1h2_aa_vs_h1h2_eur' in filename:
        return 'H1H2 AA vs H1H2 EUR'
    return filename.split('_pathway')[0]

def get_comprehensive_results(tissue_name):
    """Read the comprehensive results file for a tissue to get pathway ordering."""
    # Map tissue names to their comprehensive results files
    comp_files = {
        'astrocytes': 'gsea_comprehensive_results_astro_no_na_padj_flipped_last2.csv',
        'caudate': 'gsea_comprehensive_results_Brain_caudate_no_na_padj_flipped_last2.csv',
        'mfg': 'gsea_comprehensive_results_mfg_added_2_samples_no_na_padj_flipped_last2.csv',
        'neurons': 'gsea_comprehensive_results_neurons_no_na_padj_flipped_last2.csv'
    }
    
    try:
        comp_file = comp_files[tissue_name]
        df = pd.read_csv(comp_file)
        return df[['Pathway', 'Effect_Type']]
    except Exception as e:
        print(f"Error reading comprehensive results for {tissue_name}: {str(e)}")
        return None

def get_ordered_pathways(comp_results):
    """Order pathways based on effect types."""
    effect_type_order = [
        'up_807SNP_effect',
        'down_807SNP_effect',
        'up_haplotype_effect',
        'down_haplotype_effect',
        'up_H2_in_AFR_background',
        'down_H2_in_AFR_background'
    ]
    
    ordered_pathways = []
    for effect in effect_type_order:
        effect_pathways = comp_results[comp_results['Effect_Type'] == effect]['Pathway'].tolist()
        ordered_pathways.extend(effect_pathways)
    
    return ordered_pathways

def process_tissue_data(file_path, tissue_name):
    # Read the CSV
    df = pd.read_csv(file_path)
    
    # Clean pathway names
    df['Pathway'] = df['Pathway'].apply(clean_pathway_name)
    
    # Get comprehensive results for pathway ordering
    comp_results = get_comprehensive_results(tissue_name)
    if comp_results is not None:
        ordered_pathways = get_ordered_pathways(comp_results)
        # Clean the ordered pathway names to match
        ordered_pathways = [clean_pathway_name(p) for p in ordered_pathways]
    else:
        ordered_pathways = df['Pathway'].unique()
    
    # Extract and standardize the comparison type
    df['comparison'] = df['Source_File'].apply(get_comparison_type)
    
    # Define the desired order of comparisons
    comparison_order = [
        'AA H1H1A vs AA H1H1',
        'EUR H1H2 vs EUR H1H1',
        'AA H1H2 vs AA H1H1',
        'H1H1 AA vs H1H1 EUR',
        'H1H2 AA vs H1H2 EUR'
    ]
    
    # Create pivot table for NES values
    pivot_nes = df.pivot(
        index='Pathway',
        columns='comparison',
        values='NES'
    )
    # Reorder pathways and columns
    pivot_nes = pivot_nes.reindex(ordered_pathways)
    pivot_nes = pivot_nes.reindex(columns=comparison_order)
    
    # Create pivot table for FDR values
    pivot_fdr = df.pivot(
        index='Pathway',
        columns='comparison',
        values='FDR q-value'
    )
    # Reorder pathways and columns
    pivot_fdr = pivot_fdr.reindex(ordered_pathways)
    pivot_fdr = pivot_fdr.reindex(columns=comparison_order)
    
    return pivot_nes, pivot_fdr

def create_tissue_heatmap(tissue_file, tissue_name, transpose=False):
    """
    Create a heatmap for tissue data.
    If transpose=True, will put pathways on x-axis and comparisons on y-axis.
    """
    # Get the pivoted data
    nes_matrix, fdr_matrix = process_tissue_data(tissue_file, tissue_name.lower())
    
    # If transpose is requested, transpose the matrices
    if transpose:
        nes_matrix = nes_matrix.T
        fdr_matrix = fdr_matrix.T
    
    # Create significance annotation matrix
    stars = fdr_matrix.map(lambda x: 
        '***' if x <= 0.01 else 
        '**' if x <= 0.05 else 
        '*' if x <= 0.25 else '')
    
    # Calculate figure size to make square cells
    num_rows = len(nes_matrix)
    num_cols = len(nes_matrix.columns)
    #cell_size = 1.2  # Larger cell size
    cell_width = 1.5   # Wider cells
    cell_height = 0.6  # Shorter cells
    
    # Adjust figure size based on orientation
 #   if transpose:
        # When transposed, we need more space for pathway labels on x-axis
        # Make the figure significantly wider to accommodate the pathway names on x-axis
 #       fig_width = num_cols * cell_size + 10  # Substantially more space for pathway names
 #       fig_height = num_rows * cell_size + 5  # More space for title and y-axis labels
 #   else:
 #       # Original orientation
 #       fig_width = num_cols * cell_size + 4  # Space for labels and colorbar
 #       fig_height = num_rows * cell_size + 2  # Space for title and labels
        
    if transpose:
      fig_width = num_cols * cell_width + 10
      fig_height = num_rows * cell_height + 5
    else:
      fig_width = num_cols * cell_width + 4
      fig_height = num_rows * cell_height + 2

    
    # Create the figure with calculated dimensions
    plt.figure(figsize=(fig_width, fig_height))
    
    # Create heatmap with red-blue diverging colormap and square cells
    ax = sns.heatmap(nes_matrix,
                    cmap='RdBu_r',
                    center=0,
                    annot=False,
                    cbar_kws={'label': 'NES'},
                    robust=True)
                  #  square=True)  # Force square cells
    
    ax.set_aspect(0.4)  # Lower value = wider rectangles, adjust as needed
    
    # For transposed version, handle x-axis labels differently
    if transpose:
        # Much larger font size for transposed version
        x_fontsize = 16  # Increased from 12
        y_fontsize = 18  # Increased from 14
        
        # Pathway names are on x-axis - rotate and make them readable
        plt.xticks(rotation=45, ha='right')
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=x_fontsize, fontweight='bold')
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=y_fontsize, fontweight='bold')
        
        # Set appropriate title
        plt.title(f'{tissue_name} Pathways (Transposed)', fontsize=22)  # Increased from 18
        
        # Add extra bottom and right margin for the rotated labels
        plt.subplots_adjust(bottom=0.4, right=0.85)
        
        # Add extra space for pathway names
        plt.tight_layout(rect=[0, 0.15, 0.85, 0.95])
    else:
        # Original orientation - pathways on y-axis
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=14, rotation=45, ha='right')
        plt.title(f'{tissue_name} Pathways', fontsize=18)
    
    # Get and customize the colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)  # Increase colorbar tick font size
    cbar.set_label('NES', fontsize=16)  # Increase colorbar label font size
    
    # Add significance stars
    for i in range(len(nes_matrix.index)):
        for j in range(len(nes_matrix.columns)):
            if not pd.isna(stars.iloc[i, j]):
                # Make stars larger in transposed version
                star_fontsize = 14 if transpose else 12
                plt.text(j + 0.5, i + 0.5, stars.iloc[i, j],
                        ha='center', va='center', color='black', fontsize=star_fontsize)
    
    # Make sure the plot layout is properly adjusted but don't override the specific 
    # adjustments for transposed plots
    if not transpose:
        plt.tight_layout()
    
    # Save with higher resolution and appropriate filename
    #suffix = "_transposed" if transpose else ""
    #plt.savefig(f'{tissue_name.lower()}_heatmap{suffix}.png', dpi=300, bbox_inches='tight')
    
    
    suffix = "_transposed" if transpose else ""
    sample_suffix = "_2_added_samples" if tissue_name.lower() == 'mfg' else ""
    plt.savefig(f'{tissue_name.lower()}_heatmap{suffix}{sample_suffix}.png', dpi=300, bbox_inches='tight')

    plt.close()

# Map tissue names to their file names
tissue_files = {
    'Astrocytes': 'astro_filtered_pathways_sig_genes_from_comprehensive_gsea.csv',
    'Caudate': 'caudate_filtered_pathways_sig_genes_from_comprehensive_gsea.csv',
    'MFG': 'mfg_filtered_pathways_sig_genes_from_comprehensive_gsea_added_2_mfg_samples.csv',
    'Neurons': 'neurons_filtered_pathways_sig_genes_from_comprehensive_gsea.csv'
}

# Create both original and transposed heatmaps for each tissue
for tissue_name, file_path in tissue_files.items():
    try:
        # Create original heatmap
        create_tissue_heatmap(file_path, tissue_name, transpose=False)
        print(f"Created original heatmap for {tissue_name}")
        
        # Create transposed heatmap
        create_tissue_heatmap(file_path, tissue_name, transpose=True)
        print(f"Created transposed heatmap for {tissue_name}")
    except Exception as e:
        print(f"Error processing {tissue_name}: {str(e)}")

# Create legend figure
plt.figure(figsize=(6, 2))
plt.text(0.5, 0.5, 
         '*** FDR ≤ 0.01\n** FDR ≤ 0.05\n* FDR ≤ 0.25', 
         ha='center', va='center')
plt.axis('off')
plt.savefig('legend.png', dpi=300, bbox_inches='tight')
plt.close()

print("Done! Created all heatmaps and legend.")

def get_grouped_effect_pathways(comprehensive_files):
    """
    Get pathways grouped by effect type across all tissues.
    Returns a dictionary with three groups of tissue-specific pathways.
    """
    # Define the three effect groups
    effect_groups = {
        'snp_effects': ['up_807SNP_effect', 'down_807SNP_effect'],
        'haplotype_effects': ['up_haplotype_effect', 'down_haplotype_effect'],
        'afr_effects': ['up_H2_in_AFR_background', 'down_H2_in_AFR_background']
    }
    
    # Tissue order
    tissue_order = ['mfg', 'caudate', 'neurons', 'astro']
    
    # Dictionary to store pathways by group and tissue
    grouped_pathways = {
        'snp_effects': [],
        'haplotype_effects': [],
        'afr_effects': []
    }
    
    # Process each tissue
    for tissue_name in tissue_order:
        try:
            # Read comprehensive results file
            file_path = comprehensive_files.get(tissue_name)
            if not file_path:
                print(f"Warning: No comprehensive file found for {tissue_name}")
                continue
                
            df = pd.read_csv(file_path)
            
            # Process each effect group
            for group_name, effect_types in effect_groups.items():
                for effect_type in effect_types:
                    # Get pathways for this tissue and effect type
                    effect_pathways = df[df['Effect_Type'] == effect_type]['Pathway'].tolist()
                    
                    # Clean pathway names
                    effect_pathways = [clean_pathway_name(p) for p in effect_pathways]
                    
                    # Add to grouped pathways with tissue prefix
                    tissue_display = {
                        'mfg': 'MFG', 
                        'caudate': 'Caudate', 
                        'neurons': 'Neurons', 
                        'astro': 'Astrocytes'
                    }.get(tissue_name, tissue_name.upper())
                    
                    # Add effect type prefix
                    effect_display = {
                        'up_807SNP_effect': 'Up SNP', 
                        'down_807SNP_effect': 'Down SNP',
                        'up_haplotype_effect': 'Up Haplotype',
                        'down_haplotype_effect': 'Down Haplotype',
                        'up_H2_in_AFR_background': 'Up H2 in AFR',
                        'down_H2_in_AFR_background': 'Down H2 in AFR'
                    }.get(effect_type, effect_type)
                    
                    # Store the pathways with metadata
                    if effect_pathways:
                        grouped_pathways[group_name].append({
                            'tissue': tissue_display,
                            'effect_type': effect_display,
                            'pathways': effect_pathways
                        })
        
        except Exception as e:
            print(f"Error processing comprehensive file for {tissue_name}: {str(e)}")
    
    return grouped_pathways


def create_effect_type_heatmap(tissue_files, comprehensive_files, effect_group):
    """
    Create a combined heatmap for a specific effect group, organized by tissue.
    
    Parameters:
    -----------
    tissue_files: dict
        Dictionary mapping tissue names to their data files
    comprehensive_files: dict
        Dictionary mapping tissue names to their comprehensive results files
    effect_group: str
        One of 'snp_effects', 'haplotype_effects', or 'afr_effects'
    """
    # Get grouped pathways
    grouped_pathways = get_grouped_effect_pathways(comprehensive_files)
    
    # Get data specific to this effect group
    group_data = grouped_pathways.get(effect_group, [])
    
    if not group_data:
        print(f"No data available for effect group: {effect_group}")
        return
    
    # Determine the group title
    group_titles = {
        'snp_effects': 'SNP Effects',
        'haplotype_effects': 'Haplotype Effects',
        'afr_effects': 'H2 in AFR Background Effects'
    }
    group_title = group_titles.get(effect_group, effect_group)
    
    # Dictionary to store NES and FDR values
    all_nes_data = {}
    all_fdr_data = {}
    
    # Collect all unique pathways across all tissues for this effect group
    all_pathways = []
    for item in group_data:
        all_pathways.extend(item['pathways'])
    all_pathways = list(dict.fromkeys(all_pathways))  # Remove duplicates while preserving order
    
    # Process each tissue file
    for tissue_name, file_path in tissue_files.items():
        # Process the tissue data
        nes_matrix, fdr_matrix = process_tissue_data(file_path, tissue_name.lower())
        
        # Check which tissue+effect groups this tissue belongs to
        for item in group_data:
            # Skip if pathways list is empty
            if not item['pathways']:
                continue
                
            # Get standardized tissue name
            std_tissue = {
                'Astrocytes': 'astro',
                'Caudate': 'caudate',
                'MFG': 'mfg',
                'Neurons': 'neurons'
            }.get(tissue_name, tissue_name.lower())
            
            # Skip if this doesn't match the current tissue
            if std_tissue not in item['tissue'].lower():
                continue
            
            # For each comparison in this tissue
            for comparison in nes_matrix.columns:
                # Create a row key combining tissue, effect type, and comparison
                row_key = f"{item['tissue']} ({item['effect_type']}): {comparison}"
                
                # Create a Series indexed by all pathways, initialized with NaN
                nes_values = pd.Series(index=all_pathways, dtype='float64')
                fdr_values = pd.Series(index=all_pathways, dtype='float64')
                
                # Fill in values for pathways that exist in this tissue
                for pathway in item['pathways']:
                    if pathway in nes_matrix.index:
                        nes_values[pathway] = nes_matrix.loc[pathway, comparison]
                        fdr_values[pathway] = fdr_matrix.loc[pathway, comparison]
                
                # Store the series
                all_nes_data[row_key] = nes_values
                all_fdr_data[row_key] = fdr_values
    
    # Convert to DataFrames
    combined_nes = pd.DataFrame(all_nes_data).T
    combined_fdr = pd.DataFrame(all_fdr_data).T
    
    # Handle missing values for visualization
    combined_nes = combined_nes.fillna(0)  # Replace NaN with 0 for visualization
    
    # Create significance annotation matrix
    stars = combined_fdr.map(lambda x: 
        '***' if x <= 0.01 else 
        '**' if x <= 0.05 else 
        '*' if x <= 0.25 else '')
    stars = stars.fillna('')  # Replace NaN with empty string
    
    # Calculate figure size
    num_rows = len(combined_nes)
    num_cols = len(combined_nes.columns)
    
    # Use smaller cells for many rows/columns, larger for fewer
    #cell_size = max(0.4, min(1.0, 30 / max(num_rows, num_cols)))
    
    # Calculate dimensions, ensuring minimum reasonable size
    #fig_width = max(20, num_cols * cell_size + 12)  # Minimum 20 inches or calculated width
    #fig_height = max(15, num_rows * cell_size + 6)  # Minimum 15 inches or calculated height
    
    cell_width = 1.5
    cell_height = 0.6
    fig_width = max(20, num_cols * cell_width + 12)
    fig_height = max(8, num_rows * cell_height + 6)  # Reduced from 15

    # Create the figure
    plt.figure(figsize=(fig_width, fig_height))
    
    # Create heatmap
    ax = sns.heatmap(combined_nes,
                    cmap='RdBu_r',
                    center=0,
                    annot=False,
                    cbar_kws={'label': 'NES'},
                    robust=True)
                 #   square=True)
    
    # Format axis labels
    plt.xticks(rotation=45, ha='right')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=14, fontweight='bold')
    
    # Set title
    plt.title(f'Combined Tissues: {group_title}', fontsize=22)
    
    # Customize colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('NES', fontsize=16)
    
    # Add significance stars
    for i in range(len(combined_nes.index)):
        for j in range(len(combined_nes.columns)):
            try:
                if stars.iloc[i, j] and not pd.isna(combined_nes.iloc[i, j]) and combined_nes.iloc[i, j] != 0:
                    plt.text(j + 0.5, i + 0.5, stars.iloc[i, j],
                            ha='center', va='center', color='black', fontsize=10)
            except Exception:
                pass  # Skip if any issues with adding stars
    
    # Add extra spacing for labels
    plt.subplots_adjust(bottom=0.3, left=0.3, right=0.85)
    
    # Save the figure
    #plt.savefig(f'combined_{effect_group}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'combined_{effect_group}_heatmap_2_added_samples.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created combined heatmap for {group_title}")
    return combined_nes, combined_fdr


# Define the comprehensive results files
comprehensive_files = {
    'neurons': 'gsea_comprehensive_results_neurons_no_na_padj_flipped_last2.csv',
    'mfg': 'gsea_comprehensive_results_mfg_added_2_samples_no_na_padj_flipped_last2.csv',
    'astro': 'gsea_comprehensive_results_astro_no_na_padj_flipped_last2.csv',
    'caudate': 'gsea_comprehensive_results_Brain_caudate_no_na_padj_flipped_last2.csv'
}

# Add this after the individual tissue processing
print("\nCreating combined effect-type heatmaps...")

# Create the three effect-type heatmaps
try:
    # SNP Effects
    create_effect_type_heatmap(tissue_files, comprehensive_files, 'snp_effects')
    
    # Haplotype Effects
    create_effect_type_heatmap(tissue_files, comprehensive_files, 'haplotype_effects')
    
    # H2 in AFR Background Effects
    create_effect_type_heatmap(tissue_files, comprehensive_files, 'afr_effects')
    
    print("Created all combined effect-type heatmaps successfully")
except Exception as e:
    print(f"Error creating combined effect-type heatmaps: {str(e)}")

def create_compact_combined_heatmap(tissue_files, comprehensive_files, effect_group):
    """
    Create a combined heatmap for a specific effect group with only 5 rows (one per comparison).
    Pathways are strictly ordered by tissue and effect, preserving duplicates across tissues.
    
    Parameters:
    -----------
    tissue_files: dict
        Dictionary mapping tissue names to their data files
    comprehensive_files: dict
        Dictionary mapping tissue names to their comprehensive results files
    effect_group: str
        One of 'snp_effects', 'haplotype_effects', or 'afr_effects'
    """
    # Get data specific to this effect group
    effect_types = {
        'snp_effects': ['up_807SNP_effect', 'down_807SNP_effect'],
        'haplotype_effects': ['up_haplotype_effect', 'down_haplotype_effect'],
        'afr_effects': ['up_H2_in_AFR_background', 'down_H2_in_AFR_background']
    }.get(effect_group, [])
    
    # Determine the group title
    group_titles = {
        'snp_effects': 'SNP Effects',
        'haplotype_effects': 'Haplotype Effects',
        'afr_effects': 'H2 in AFR Background Effects'
    }
    group_title = group_titles.get(effect_group, effect_group)
    
    # Dictionary to store NES and FDR values for each comparison
    all_nes_data = {}
    all_fdr_data = {}
    
    # Define tissue order and effect type order - explicit ordering
    tissue_order = ['mfg', 'caudate', 'neurons', 'astro']
    effect_prefix_map = {
        'up_807SNP_effect': 'Up',
        'down_807SNP_effect': 'Down',
        'up_haplotype_effect': 'Up',
        'down_haplotype_effect': 'Down',
        'up_H2_in_AFR_background': 'Up',
        'down_H2_in_AFR_background': 'Down'
    }
    
    # Create a list to store pathways in the exact tissue+effect order,
    # including duplicates across different tissues
    ordered_pathways_with_tissue = []
    
    # Process each tissue to get its pathways in order
    for tissue_abbr in tissue_order:
        file_path = comprehensive_files.get(tissue_abbr)
        if not file_path:
            print(f"Warning: No comprehensive file found for {tissue_abbr}")
            continue
            
        try:
            # Read the file
            df = pd.read_csv(file_path)
            
            # Get tissue display name
            tissue_display = {
                'mfg': 'MFG',
                'caudate': 'Caudate',
                'neurons': 'Neurons',
                'astro': 'Astrocytes'
            }.get(tissue_abbr, tissue_abbr.upper())
            
            # Process effect types in order (Up then Down)
            for effect_type in effect_types:
                # Skip if this effect type doesn't apply to this effect group
                if effect_type not in effect_prefix_map:
                    continue
                    
                # Get effect display name
                effect_display = effect_prefix_map[effect_type]
                
                # Get pathways for this tissue and effect type
                pathways_df = df[df['Effect_Type'] == effect_type]
                
                # Process each pathway
                for _, row in pathways_df.iterrows():
                    pathway = clean_pathway_name(row['Pathway'])
                    
                    # Create a unique identifier for this tissue+pathway
                    tissue_pathway = f"{tissue_display}:{pathway}"
                    
                    # Add to ordered list with metadata
                    ordered_pathways_with_tissue.append({
                        'tissue': tissue_display,
                        'tissue_abbr': tissue_abbr,
                        'effect': effect_display,
                        'effect_type': effect_type,
                        'pathway': pathway,
                        'tissue_pathway': tissue_pathway  # Unique identifier
                    })
        
        except Exception as e:
            print(f"Error processing comprehensive file for {tissue_abbr}: {str(e)}")
    
    # Get unique tissue_pathway identifiers to use as columns
    tissue_pathway_cols = [item['tissue_pathway'] for item in ordered_pathways_with_tissue]
    pathway_cols = [item['pathway'] for item in ordered_pathways_with_tissue]
    
    # Define the desired order of comparisons
    comparison_order = [
        'AA H1H1A vs AA H1H1',
        'EUR H1H2 vs EUR H1H1',
        'AA H1H2 vs AA H1H1',
        'H1H1 AA vs H1H1 EUR',
        'H1H2 AA vs H1H2 EUR'
    ]
    
    # Initialize DataFrames with comparisons as rows and tissue_pathways as columns
    nes_df = pd.DataFrame(index=comparison_order, columns=tissue_pathway_cols)
    fdr_df = pd.DataFrame(index=comparison_order, columns=tissue_pathway_cols)
    
    # Fill in the data from each tissue file
    for tissue_name, file_path in tissue_files.items():
        # Get standard tissue abbreviation
        std_tissue = {
            'Astrocytes': 'astro',
            'Caudate': 'caudate',
            'MFG': 'mfg',
            'Neurons': 'neurons'
        }.get(tissue_name, tissue_name.lower())
        
        # Process the tissue file
        nes_matrix, fdr_matrix = process_tissue_data(file_path, tissue_name.lower())
        
        # For each pathway in our ordered list that belongs to this tissue
        for item in ordered_pathways_with_tissue:
            if item['tissue_abbr'] != std_tissue:
                continue
                
            pathway = item['pathway']
            tissue_pathway = item['tissue_pathway']
            
            # For each comparison
            for comparison in comparison_order:
                if comparison in nes_matrix.columns and pathway in nes_matrix.index:
                    nes_df.loc[comparison, tissue_pathway] = nes_matrix.loc[pathway, comparison]
                    fdr_df.loc[comparison, tissue_pathway] = fdr_matrix.loc[pathway, comparison]
    
    # Replace NaN with 0 for visualization
    nes_df = nes_df.fillna(0)
    
    # Create significance annotation matrix
    stars_df = fdr_df.map(lambda x: 
        '***' if x <= 0.01 else 
        '**' if x <= 0.05 else 
        '*' if x <= 0.25 else '')
    stars_df = stars_df.fillna('')
    
    # Calculate figure size
    num_rows = len(nes_df)  # Should be 5 (one per comparison)
    num_cols = len(nes_df.columns)  # Number of tissue-specific pathways
    
    # Calculate cell size - smaller for many columns, larger for fewer
   # cell_size = max(0.2, min(0.8, 30 / num_cols))
    
    # Calculate dimensions with minimum reasonable size
   # fig_width = max(24, num_cols * cell_size + 8)
   # fig_height = max(10, num_rows * cell_size + 8)
    
    
     # Remove square=True and adjust dimensions similarly
    cell_width = 1.5
    cell_height = 0.6
    fig_width = max(24, num_cols * cell_width + 8)
    fig_height = max(6, num_rows * cell_height + 4)  # Reduced height

    
    # Create the figure
    plt.figure(figsize=(fig_width, fig_height))
    
    # Create heatmap
    ax = sns.heatmap(nes_df,
                    cmap='RdBu_r',
                    center=0,
                    annot=False,
                    cbar_kws={'label': 'NES'},
                    robust=True)
                #    square=True)
    
    # Replace tissue_pathway IDs with just pathway names for display
    ax.set_xticklabels(pathway_cols, rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=14, fontweight='bold')
    
    # Set title
    plt.title(f'Combined Tissues: {group_title} (Compact View)', fontsize=22)
    
    # Customize colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('NES', fontsize=16)
    
    # Add significance stars
    for i in range(len(nes_df.index)):
        for j in range(len(nes_df.columns)):
            try:
                if stars_df.iloc[i, j] and nes_df.iloc[i, j] != 0:
                    plt.text(j + 0.5, i + 0.5, stars_df.iloc[i, j],
                            ha='center', va='center', color='black', fontsize=8)
            except Exception:
                pass  # Skip if any issues with adding stars
    
    # Add extra spacing for labels with slanted text
    plt.subplots_adjust(bottom=0.3, right=0.85)
    
    # Save the figure
   # plt.savefig(f'compact_combined_{effect_group}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'compact_combined_{effect_group}_heatmap_2_added_samples.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created compact combined heatmap for {group_title} preserving duplicates across tissues")
    print(f"Generated with {len(tissue_pathway_cols)} columns including duplicate pathways")
    return nes_df, fdr_df

# Add this after the individual tissue processing
print("\nCreating compact combined effect-type heatmaps with duplicates preserved...")

try:
    # SNP Effects (Compact)
    create_compact_combined_heatmap(tissue_files, comprehensive_files, 'snp_effects')
    
    # Haplotype Effects (Compact)
    create_compact_combined_heatmap(tissue_files, comprehensive_files, 'haplotype_effects')
    
    # H2 in AFR Background Effects (Compact)
    create_compact_combined_heatmap(tissue_files, comprehensive_files, 'afr_effects')
    
    print("Created all compact combined effect-type heatmaps successfully")
except Exception as e:
    print(f"Error creating compact combined effect-type heatmaps: {str(e)}")