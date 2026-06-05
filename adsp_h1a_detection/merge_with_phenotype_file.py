#!/usr/bin/env python3
import pandas as pd

# ===============================================================
# CONFIG — edit these paths for your environment.
# ===============================================================
# Phenotype table (tab-separated). Must contain `sample.id` plus the
# phenotype/covariate columns used downstream (e.g. diagnosis, sex, apoe,
# sequencing_center, sequencing_platform).
PHENOTYPE_TSV = "path/to/adsp_phenotypes.tsv"

# Global-ancestry / admixture table (tab-separated). Must contain `IID`
# plus the ancestry-fraction columns and `admixture_super_pop_max`.
ANCESTRY_TSV = "path/to/admixture_global_ancestry.tsv"

# Sample-annotation table with the PCA covariate columns V1..Vn
# (tab-separated, keyed on `sample.id`).
COVARIATE_TSV = "path/to/adsp_sample_annot_merged_covars.tsv"
# ===============================================================

def merge_phenotype_and_ancestry_data():
    """
    Complete pipeline to merge genotype data with phenotype and ancestry information.
    """
    # File paths
    genotype_file = 'results/mapt_genotypes_v10_with_agreement.csv'
    phenotype_file = PHENOTYPE_TSV
    ancestry_file = ANCESTRY_TSV

    # -----------------------------
    # Step 1: Load and prepare genotype data
    # -----------------------------
    print("Loading genotype data...")
    df = pd.read_csv(genotype_file, low_memory=False)
    
    # Ensure sample_name_cut exists
    if 'sample_name_cut' not in df.columns:
        df['sample_name_cut'] = df['Sample'].apply(lambda x: "-".join(x.split('-')[:3]))
    
    # -----------------------------
    # Step 2: Merge with phenotype data
    # -----------------------------
    print("Loading and merging phenotype data...")
    pheno_df = pd.read_csv(phenotype_file, sep='\t', low_memory=False)
    merged_df = pd.merge(df, pheno_df, left_on='Sample', right_on='sample.id', how='left')
    
    # Print phenotype merge statistics
    total_samples = len(merged_df)
    num_with_pheno = merged_df['sample.id'].notnull().sum()
    print(f"\nPhenotype Merge Statistics:")
    print(f"Total samples in genotype file: {total_samples}")
    print(f"Number of samples with phenotype metadata: {num_with_pheno}")
    
    # Agreement group statistics
    agree_no_pheno = merged_df[(merged_df['Agreement'] == 'Agree') & (merged_df['sample.id'].isna())].shape[0]
    disagree_no_pheno = merged_df[(merged_df['Agreement'] == 'Disagree') & (merged_df['sample.id'].isna())].shape[0]
    print(f"Samples with 'Agree' lacking phenotype data: {agree_no_pheno}")
    print(f"Samples with 'Disagree' lacking phenotype data: {disagree_no_pheno}")
    
    # -----------------------------
    # Step 3: Add ancestry data
    # -----------------------------
    print("\nLoading and merging ancestry data...")
    ancestry_df = pd.read_csv(ancestry_file, sep='\t', low_memory=False)
    
    # Select only needed ancestry columns
    ancestry_columns = ['IID', 'AFR', 'AMR', 'EAS', 'EUR', 'OCE', 'SAS','admixture_super_pop_max', 'admixture_cluster_max']
    ancestry_subset = ancestry_df[ancestry_columns]
    
    # Merge with ancestry data
    final_df = pd.merge(
        merged_df,
        ancestry_subset,
        left_on='Sample',
        right_on='IID',
        how='left'
    )
    
    # Print ancestry merge statistics
    samples_with_ancestry = final_df[['AFR', 'AMR', 'EAS', 'EUR', 'OCE', 'SAS','admixture_super_pop_max', 'admixture_cluster_max']].notna().any(axis=1).sum()
    print(f"\nAncestry Merge Statistics:")
    print(f"Total samples: {total_samples}")
    print(f"Samples with ancestry data: {samples_with_ancestry}")
    print(f"Samples without ancestry data: {total_samples - samples_with_ancestry}")
    print(f"Percentage of samples with ancestry data: {(samples_with_ancestry/total_samples)*100:.2f}%")
    
    # -----------------------------
    # Step 4: Save final merged dataset
    # -----------------------------
    output_file = 'results/mapt_genotypes_v10_with_agreement_and_pheno_and_global_ancestry.csv'
    print(f"\nSaving complete merged dataset to: {output_file}")
    
    # Drop the IID column since it's redundant with Sample
    final_df = final_df.drop('IID', axis=1, errors='ignore')
    
    # Save to CSV
    final_df.to_csv(output_file, index=False)
    print("Merge complete!")

if __name__ == "__main__":
    merge_phenotype_and_ancestry_data()

# -----------------------------
# Step 5: Attach the PCA covariate columns (V1..Vn) from the sample-annotation
# table to produce the final modeling input.
# -----------------------------
tsv_file = COVARIATE_TSV
csv_file = "results/mapt_genotypes_v10_with_agreement_and_pheno_and_global_ancestry.csv"
output_file = "results/mapt_genotypes_with_covar_merged.csv"

# Read the files
print(f"Reading TSV file: {tsv_file}")
tsv_data = pd.read_csv(tsv_file, sep='\t')

print(f"Reading CSV file: {csv_file}")
csv_data = pd.read_csv(csv_file)

# Extract only the sample.id and V1-V32 columns from TSV
v_columns = [col for col in tsv_data.columns if col.startswith('V')]
print(f"Found {len(v_columns)} V columns: {', '.join(v_columns)}")

# Handle the V25V26 issue if it exists
if 'V25V26' in v_columns:
    print("Found problematic column V25V26, renaming to V25")
    tsv_data = tsv_data.rename(columns={'V25V26': 'V25'})
    v_columns = [col for col in tsv_data.columns if col.startswith('V')]

# Create a subset of the TSV with only sample.id and V columns
tsv_subset = tsv_data[['sample.id'] + v_columns]

# Verify both dataframes before merging
print(f"TSV data shape: {tsv_data.shape}")
print(f"CSV data shape: {csv_data.shape}")
print(f"TSV subset shape: {tsv_subset.shape}")

# Check if sample.id exists in both files
if 'sample.id' not in tsv_subset.columns:
    print("Error: 'sample.id' column not found in TSV file")
    exit(1)
if 'sample.id' not in csv_data.columns:
    print("Error: 'sample.id' column not found in CSV file")
    exit(1)

# Merge on sample.id
print("Merging files on sample.id column...")
merged_data = pd.merge(csv_data, tsv_subset, on='sample.id', how='left')

# Check for NaN values after merging (indicates samples that didn't match)
nan_count = merged_data[v_columns].isna().sum().sum()
print(f"Number of NaN values in V columns after merge: {nan_count}")
print(f"Merged data shape: {merged_data.shape}")

# Save the merged data
print(f"Saving merged data to: {output_file}")
merged_data.to_csv(output_file, index=False)

print("Merge complete!")