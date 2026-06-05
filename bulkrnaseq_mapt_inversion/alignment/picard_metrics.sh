#!/bin/bash

# Setting environment variables
#JAVA_FLAGS="-Xmx4g -Xms64m"
REF_FILE="path/to/BulkRNASeq/References/v38/GRCh38.p13.genome.fa"

TX_FLAT_FILE="path/to/BulkRNASeq/References/v38/gencode.v38.refFlat"
# Convert GTF to genePred format
# Convert genePred to RefFlat format

# Cleanup temporary files


#conda activate gatk_env
 
#conda activate gatk_env

RIBOSOMAL_INTERVALS="path/to/BulkRNASeq/References/v38/gencode.v38.annotation_rrna.ribosomal_intervals.interval_list"
RNA_STRAND=SECOND_READ_TRANSCRIPTION_STRAND
PICARD_INPUT="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/astro/picard_params_redo.txt"
LOG_DIR="path/to/BulkRNASeq/Ancestry_analysis/no_wasp_outputs/astro/picard/logs"

# Ensure the log directory exists
mkdir -p "${LOG_DIR}"

# Read from params.txt
while read -r IN_BAM INTERVAL_LIST METRICS_FILE; do
    START=$SECONDS

    SAMPLE_ID=$(basename "$IN_BAM" _Aligned.sortedByCoord.out.bam)
    # Extract directory from METRICS_FILE for output
    OUT_DIR=$(dirname "${METRICS_FILE}")
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.log"
    #PICARD_DONE="${METRICS_FILE}.picard.done"

    echo "Processing sample: $SAMPLE_ID" | tee -a "${SAMPLE_LOG}"
    echo "Start time: $(date)" | tee -a "${SAMPLE_LOG}"

        # CollectAlignmentSummaryMetrics
        ALIGNMENT_METRICS_OUT="${OUT_DIR}/${SAMPLE_ID}.alignment_metrics.txt"
        echo "Running CollectAlignmentSummaryMetrics..." | tee -a "${SAMPLE_LOG}"
        java -jar path/to/picard.jar CollectAlignmentSummaryMetrics \
        R="$REF_FILE" \
        I="$IN_BAM" \
        O="$ALIGNMENT_METRICS_OUT" \
        EXPECTED_PAIR_ORIENTATIONS=FR \
        EXPECTED_PAIR_ORIENTATIONS=RF 2>> "${SAMPLE_LOG}"

        # CollectRnaSeqMetrics
        RNA_SEQ_METRICS_OUT="${OUT_DIR}/${SAMPLE_ID}.rna_seq_metrics.txt"
        echo "Running CollectRnaSeqMetrics..." | tee -a "${SAMPLE_LOG}"
        java -jar path/to/picard.jar CollectRnaSeqMetrics \
        I="$IN_BAM" \
        O="$RNA_SEQ_METRICS_OUT" \
        R="$REF_FILE" \
        STRAND_SPECIFICITY="$RNA_STRAND" \
        REF_FLAT="$TX_FLAT_FILE" \
        RIBOSOMAL_INTERVALS="$RIBOSOMAL_INTERVALS" 2>> "${SAMPLE_LOG}"

        # Collect GC Bias Metrics
        GC_METRICS_OUT="${OUT_DIR}/${SAMPLE_ID}.gc_bias_metrics.txt"
        GC_SUMMARY_OUT="${OUT_DIR}/${SAMPLE_ID}.gc_bias_metrics.summary.txt"
        CHART_OUT="${OUT_DIR}/${SAMPLE_ID}.gc_bias.chart.pdf"
        echo "Running CollectGcBiasMetrics..." | tee -a "${SAMPLE_LOG}"
        java -jar path/to/picard.jar CollectGcBiasMetrics \
        R="$REF_FILE" \
        I="$IN_BAM" \
        O="$GC_METRICS_OUT" \
        CHART_OUTPUT="$CHART_OUT" \
        SUMMARY_OUTPUT="$GC_SUMMARY_OUT" 2>> "${SAMPLE_LOG}"

        # Mark Duplicates
        DUP_METRICS_OUT="${OUT_DIR}/${SAMPLE_ID}.duplication_metrics.txt"
        DUP_BAM_OUT="${OUT_DIR}/${SAMPLE_ID}.marked_dup.bam"
        echo "Running MarkDuplicates..." | tee -a "${SAMPLE_LOG}"
        java -jar path/to/picard.jar MarkDuplicates \
        I="$IN_BAM" \
        O="$DUP_BAM_OUT" \
        M="$DUP_METRICS_OUT" 2>> "${SAMPLE_LOG}"

    END=$SECONDS
    echo "Time taken: $((END - START)) seconds" | tee -a "${SAMPLE_LOG}"
    echo "--------------------------" | tee -a "${SAMPLE_LOG}"

done < "$PICARD_INPUT"