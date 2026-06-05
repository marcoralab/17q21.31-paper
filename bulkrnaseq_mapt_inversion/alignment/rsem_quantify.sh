#!/bin/bash

# Define directories and parameters
LOG_DIR="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/astro/rsem/logs"
RSEM_IDX="path/to/BulkRNASeq/References/v38/rsem_reference_v38"
RSEM_INPUT_FILE="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/astro/rsem_params_v38.txt"
RSEM_FMT_BIN="convert-sam-for-rsem"
RSEM_QUANT_BIN="rsem-calculate-expression"
SAMTOOLS_BIN="samtools"
THREADS=4
STRANDEDNESS=none

# Ensure the log directory exists
mkdir -p "${LOG_DIR}"

# Read from the parameter file
while read -r SAMPLE_ID IN_TX_BAMS OUT_SAMPLE_TX_BAM OUT_FILE; do
    START=$SECONDS  # Initialize START variable for each sample processing

    # Setup variables
    OUT_DIR=$(dirname "${OUT_SAMPLE_TX_BAM}")
    mkdir -p "${OUT_DIR}"
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.log"
    echo "Processing sample: $SAMPLE_ID" >> "${SAMPLE_LOG}"

    OUT_PFX=$(echo ${OUT_SAMPLE_TX_BAM} | sed 's/\.tx\.bam//')
    OUT_FILE_DONE="$OUT_FILE.done"

    if [ ! -f $OUT_FILE_DONE ]; then
        if [ -f $IN_TX_BAMS ]; then
            echo "Copying BAM file for $SAMPLE_ID" >> "${SAMPLE_LOG}"
            cp $IN_TX_BAMS $OUT_SAMPLE_TX_BAM
            copy_success=$?

            # Check if the BAM file was created
            if [ ! -f "$OUT_SAMPLE_TX_BAM" ]; then
                echo "ERROR: BAM file not created for $SAMPLE_ID" | tee -a "${SAMPLE_LOG}"
                exit 1
            fi
        else
            echo "Input BAM file $IN_TX_BAMS does not exist, check ${SAMPLE_LOG}" >> "${SAMPLE_LOG}"
            exit 1
        fi

        echo "Calculating expression for $SAMPLE_ID" >> "${SAMPLE_LOG}"
        $RSEM_QUANT_BIN --paired-end \
        --strandedness $STRANDEDNESS \
        --num-threads $THREADS \
        --alignments \
        --estimate-rspd \
        --no-bam-output \
        --seed 20161010 \
        $OUT_SAMPLE_TX_BAM \
        $RSEM_IDX \
        $OUT_PFX 2>>"${SAMPLE_LOG}"

        quant_success=$?

        if [ $copy_success -eq 0 ] && [ $quant_success -eq 0 ]; then
            rm $OUT_SAMPLE_TX_BAM
            touch $OUT_FILE_DONE
        else
            echo "Error in processing $SAMPLE_ID, check logs." >> "${SAMPLE_LOG}"
            exit 1
        fi
    else
        echo "$OUT_FILE_DONE exists, skipping $SAMPLE_ID." >> "${SAMPLE_LOG}"
    fi

    # If this was faster than 10 minutes, sleep for the rest.
    END=$SECONDS
    ELAPSED=$((END-START))
    if [ "$ELAPSED" -lt 600 ]; then
        TOSLEEP=$((600 - ELAPSED))
        sleep $TOSLEEP
    fi

done < "$RSEM_INPUT_FILE"