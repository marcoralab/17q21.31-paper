#!/bin/bash

# Directory containing FASTQ files
FASTQ_DIR="path/to/fastq_dir/" 

# Output directory for STAR
OUT_DIR="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/microglia"

# Base name for read group IDs, change this according to your needs
RG_BASE="RG"

# File to store parameters
PARAMS_FILE="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/microglia/star_params.txt"
> "$PARAMS_FILE"  # Clear the file to start fresh

# Loop through FASTQ files and create parameters
for R1_FILE in ${FASTQ_DIR}/*Microglia*_R1_001.fastq.gz; do
    # Extract basic sample name by removing path and file extension
    SAMPLE_NAME=$(basename "${R1_FILE}" "_R1_001.fastq.gz")
    
    # Corresponding R2 file
    R2_FILE="${FASTQ_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"
    
    # Construct the read group ID (modify this line to match your naming convention)
    IN_RG="${RG_BASE}_${SAMPLE_NAME}_L001"  # Example: RG_SampleID_L001
    
    # Construct the output path for BAM files
    OUT_BAM="${OUT_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"
    
    # Write to params file
    echo "${SAMPLE_NAME} ${R1_FILE} ${R2_FILE} ${IN_RG} ${OUT_BAM}" >> "$PARAMS_FILE"
done

echo "Parameter file created: $PARAMS_FILE"