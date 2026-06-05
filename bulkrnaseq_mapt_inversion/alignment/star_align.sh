#!/bin/bash
# Define directories and parameters
QMODE="TranscriptomeSAM"
MULT_ORDER="Random"
OUT_MULT_MAP_READS=1
VCF="path/to/BulkRNASeq/Corces_2020_NatGenet_2366H1vsH2SNPs.vcf"
LOG_DIR="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/astro/logs"
STAR_INDEX="path/to/BulkRNASeq/References/v38/star_genome_index_v38"
STAR_INPUT_FILE="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/astro/star_params.txt"

# Ensure the log directory exists
mkdir -p "${LOG_DIR}"
 
# Read from the parameter file
while read -r SAMPLE_ID IN_FQ1 IN_FQ2 IN_RG OUT_BAM; do
    # Extract directory from OUT_BAM for output
    OUT_DIR=$(dirname "${OUT_BAM}")
    # Create directory if it doesn't exist
    mkdir -p "${OUT_DIR}"
    # Setup log file path
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.log"
    
    # Logging
    echo "Processing ${SAMPLE_ID}" | tee -a "${SAMPLE_LOG}"
    # Construct read group line
    RG_LINE="ID:${IN_RG} SM:${SAMPLE_ID} PL:ILLUMINA"
    # Run STAR alignment
    echo "Running STAR Alignment for ${SAMPLE_ID}" | tee -a "${SAMPLE_LOG}"
    START_TIME=$(date +%s)
    
    STAR --runMode alignReads \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${IN_FQ1} ${IN_FQ2} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${OUT_DIR}/${SAMPLE_ID}_ \
         --runThreadN 4 \
         --outMultimapperOrder ${MULT_ORDER} \
         --outSAMmultNmax ${OUT_MULT_MAP_READS} \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within KeepPairs \
         --outBAMcompression 6 \
         --bamRemoveDuplicatesType UniqueIdentical \
         --outSJfilterReads Unique \
         --quantMode ${QMODE} \
         --outSAMattrRGline ${RG_LINE}

    # Time and check
    END_TIME=$(date +%s)
    echo "STAR Alignment for ${SAMPLE_ID} took $((END_TIME - START_TIME)) seconds" | tee -a "${SAMPLE_LOG}"
    
    # Check if STAR output was successfully created
    if [ ! -f "${OUT_BAM}" ] || [ ! -s "${OUT_BAM}" ]; then
        echo "ERROR: STAR alignment failed for ${SAMPLE_ID}. Expected output file ${OUT_BAM} not found or is empty." | tee -a "${SAMPLE_LOG}"
        exit 1
    fi
done < "$STAR_INPUT_FILE"