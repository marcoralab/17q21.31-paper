#!/usr/bin/env python3

import subprocess
import os
import logging
import pandas as pd

# ===============================================================
# CONFIG — edit for your environment.
# ===============================================================
# Bgzipped, tabix-indexed chr17 VCF from the ADSP r4 WGS release
# (controlled-access; obtain via NIAGADS/dbGaP). Requires `bcftools`
# on your PATH.
VCF_FILE = "path/to/adsp_r4_wgs_chr17.vcf.bgz"
# ===============================================================

# Setup basic logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(levelname)s - %(message)s')

# MAPT H1/H2 tag SNPs
snp_positions = {
    'rs8070723': '46003698',
    'rs1052553': '45996523'
}

# Create output directory
os.makedirs('results', exist_ok=True)

# Store genotypes for each SNP
genotypes = {}

for rs_id, pos in snp_positions.items():
    logging.info(f"Extracting genotypes for {rs_id}")
    
    # Use bcftools view to get the full VCF line
    cmd = f"bcftools view -H -r chr17:{pos} {VCF_FILE}"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        vcf_line = result.stdout.strip()
        
        # Split VCF line into fields
        fields = vcf_line.split('\t')
        
        # Get sample genotypes (starts after FORMAT field)
        format_idx = fields.index('GT')
        genotype_data = fields[format_idx + 1:]
        
        # Get sample IDs (from bcftools query)
        cmd_samples = f"bcftools query -l {VCF_FILE}"
        result_samples = subprocess.run(cmd_samples, shell=True, capture_output=True, text=True)
        
        if result_samples.returncode == 0:
            sample_ids = result_samples.stdout.strip().split('\n')
            
            # Store genotypes
            for sample_id, gt in zip(sample_ids, genotype_data):
                if gt != "./.":  # Skip missing genotypes
                    if sample_id not in genotypes:
                        genotypes[sample_id] = {}
                    genotypes[sample_id][rs_id] = gt
                    
            logging.info(f"Processed {rs_id}: found {len(genotype_data)} genotypes")
        else:
            logging.error(f"Error getting sample IDs: {result_samples.stderr}")
    else:
        logging.error(f"Error running bcftools: {result.stderr}")
        continue

# Convert to DataFrame
if genotypes:
    # Create list of records
    records = []
    for sample_id, gts in genotypes.items():
        if len(gts) == len(snp_positions):  # Only include samples with both SNPs
            records.append({
                'Sample': sample_id,
                'rs8070723': gts.get('rs8070723', './.'
),
                'rs1052553': gts.get('rs1052553', './.')
            })
    
    # Create DataFrame
    df = pd.DataFrame(records)
    
    # Show first few rows
    print("\nFirst few rows of data:")
    print(df.head())
    
    # Save to CSV
    output_file = 'results/mapt_genotypes_v10.csv'
    df.to_csv(output_file, index=False)
    
    # Print summary
    print(f"\nTotal samples: {len(df)}")
    print("\nGenotype counts:")
    for snp in snp_positions.keys():
        print(f"\n{snp}:")
        print(df[snp].value_counts())
    
    print(f"\nResults saved to: {output_file}")
else:
    logging.error("No genotype data was extracted")