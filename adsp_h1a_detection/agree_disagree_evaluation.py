#!/usr/bin/env python3
import pandas as pd

# Load the CSV file you previously generated
input_file = 'results/mapt_genotypes_v10.csv'
df = pd.read_csv(input_file)

# Add an 'Agreement' column: 'Agree' if both SNP genotypes are identical, otherwise 'Disagree'
df['Agreement'] = df.apply(
    lambda row: 'Agree' if row['rs8070723'] == row['rs1052553'] else 'Disagree',
    axis=1
)

# Optionally, add a 'Haplotype' column (here we simply join the two genotype strings)
df['Haplotype'] = df['rs8070723'] + "_" + df['rs1052553']

# Create a 'sample_name_cut' column by splitting the Sample string on '-' and joining the first three parts
df['sample_name_cut'] = df['Sample'].apply(lambda x: "-".join(x.split('-')[:3]))

# Save the updated DataFrame to a new CSV file
output_file = 'results/mapt_genotypes_v10_with_agreement.csv'
df.to_csv(output_file, index=False)

# Print summary statistics
total_samples = len(df)
agree_count = (df['Agreement'] == 'Agree').sum()
disagree_count = (df['Agreement'] == 'Disagree').sum()

print(f"Updated CSV with Agreement, Haplotype, and sample_name_cut columns saved to: {output_file}")
print(f"Total samples: {total_samples}")
print(f"Agree: {agree_count} ({agree_count/total_samples*100:.1f}%)")
print(f"Disagree: {disagree_count} ({disagree_count/total_samples*100:.1f}%)")

# If there are any disagreeing samples, print a detailed summary of their genotype combinations,
# excluding the specified combinations.
if disagree_count > 0:
    # Filter out the disagreeing samples
    disagree_df = df[df['Agreement'] == 'Disagree']
    
    # Define the genotype combinations to exclude as tuples: (rs8070723, rs1052553)
    exclude_combinations = [
        ('0/0', '0/1'),
        ('0|1', '0/1'),
        ('0/1', '1/1'),
        ('0/0', '0|0'),
        ('0/0', '0|1')
    ]
    
    # Exclude rows matching the above genotype combinations using .iloc for positional indexing.
    disagree_df_filtered = disagree_df[
        ~disagree_df[['rs8070723', 'rs1052553']].apply(
            lambda x: (x.iloc[0], x.iloc[1]), axis=1
        ).isin(exclude_combinations)
    ]
    
    # Group the remaining disagreeing samples by genotype combination and count
    disagree_summary = disagree_df_filtered.groupby(
        ['rs8070723', 'rs1052553']
    ).size().reset_index(name='Count')
    
    print("\nSummary of Disagreeing Genotype Combinations (excluding specified ones):")
    if disagree_summary.empty:
        print("No disagreeing genotype combinations remain after excluding the specified ones.")
    else:
        for idx, row in disagree_summary.iterrows():
            print(f"rs8070723: {row['rs8070723']}, rs1052553: {row['rs1052553']} -> Count: {row['Count']}")
else:
    print("No disagreeing samples found.")

