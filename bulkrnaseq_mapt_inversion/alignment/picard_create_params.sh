#!/bin/bash

# Base directory for input BAM files
IN_BAM_DIR="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/neurons" # Directory containing your input BAM files

# Directory for Picard results
PICARD_OUT_DIR="${IN_BAM_DIR}/picard"
mkdir -p "$PICARD_OUT_DIR" 

# File to store Picard parameters
PICARD_PARAMS_FILE="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/neurons/picard_params_redo.txt"
> "$PICARD_PARAMS_FILE"  # Clear the file to start fresh

# Loop through input BAM files to create parameters
for IN_BAM in ${IN_BAM_DIR}/*_Aligned.sortedByCoord.out.bam; do
    # Extract basic sample name by removing path and file extension
    SAMPLE_NAME=$(basename "${IN_BAM}" ".bam")
    
    # Define output paths
    INTERVAL_LIST="${PICARD_OUT_DIR}/${SAMPLE_NAME}.interval_list"
    METRICS_FILE="${PICARD_OUT_DIR}/${SAMPLE_NAME}_metrics.txt"

    # Write to Picard params file
    echo "${IN_BAM} ${INTERVAL_LIST} ${METRICS_FILE}" >> "$PICARD_PARAMS_FILE"
done

echo "Picard parameter file created: $PICARD_PARAMS_FILE"
