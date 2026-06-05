#!/bin/bash

# Base directory for input BAM files
IN_BAM_DIR="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/neurons"
# Output base directory
OUT_DIR_BASE="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/neurons" # Your desired output base directory

# Directory for RSEM results
#RSEM_OUT_DIR="${OUT_DIR_BASE}/results_rerun/rsem"
RSEM_OUT_DIR="${OUT_DIR_BASE}/rsem"
mkdir -p "$RSEM_OUT_DIR" 

# File to store RSEM parameters
RSEM_PARAMS_FILE="${IN_BAM_DIR}/rsem_params_v38.txt"
> "$RSEM_PARAMS_FILE"  # Clear the file to start fresh

# Loop through transcriptome BAM files to create parameters
for IN_BAM in ${IN_BAM_DIR}/*_Aligned.toTranscriptome.out.bam; do
    # Extract basic sample name by removing path and file extension
    SAMPLE_NAME=$(basename "${IN_BAM}" "_Aligned.toTranscriptome.out.bam")
    
    # Define output paths
    OUT_SAMPLE_TX_BAM="${RSEM_OUT_DIR}/${SAMPLE_NAME}.tx.bam"
    OUT_FILE="${RSEM_OUT_DIR}/${SAMPLE_NAME}.quant"

    # Write to RSEM params file
    echo "${SAMPLE_NAME} ${IN_BAM} ${OUT_SAMPLE_TX_BAM} ${OUT_FILE}" >> "$RSEM_PARAMS_FILE"
done

echo "RSEM parameter file created: $RSEM_PARAMS_FILE"
