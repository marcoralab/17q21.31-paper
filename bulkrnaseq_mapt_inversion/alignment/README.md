# Alignment & quantification (STAR → Picard → RSEM)

FASTQ → per-sample RSEM gene counts. Each stage has a parameter-file generator
and a job script that consumes it; edit the path variables at the top of each.

| Stage | Generator | Job script | Output |
|-------|-----------|-----------|--------|
| Align | `create_star_params.sh` | `star_align.sh` | `*_Aligned.sortedByCoord.out.bam`, `*_Aligned.toTranscriptome.out.bam` |
| QC    | `picard_create_params.sh` | `picard_metrics.sh` | Picard alignment / RNA-seq / GC / duplication metrics |
| Quant | `rsem_create_params.sh` | `rsem_quantify.sh` | `*.genes.results`, `*.isoforms.results` |

Run order per cohort: `create_star_params.sh` → `star_align.sh` →
`picard_create_params.sh` → `picard_metrics.sh` → `rsem_create_params.sh` →
`rsem_quantify.sh`.

The same scripts were run separately per cell type (astrocytes / neurons /
microglia) and per brain region (DLPFC / MFG / caudate) — only the FASTQ/output
directories change. Shown here is one representative set.

Requires `STAR`, Picard (`picard.jar`), `RSEM`, `samtools`, and the matching
GRCh38 (GENCODE v38) references: STAR index, RSEM reference, genome FASTA,
refFlat, and rRNA interval list.
