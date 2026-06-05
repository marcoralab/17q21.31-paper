# Iso-seq analysis: exon-level isoform calling and tissue-stratified plots

This directory contains the iso-seq code used to produce the publication
figures in:

- `tissue_plots/sample_variability_plots/sample_variability_<tissue>.png`
- `tissue_plots/read_count_over_total_sample_read_count/<tissue>/ancestry_haplotype_all_isoforms_heatmap.png`

The analysis covers wild-type (WT) samples across EUR and AA ancestry.

## Files

| File                                          | Purpose                                                                                                |
|-----------------------------------------------|--------------------------------------------------------------------------------------------------------|
| `final_process_all_isoseq_data.py`            | Master processing script. Reads per-sample BED + classification files, maps to MAPT exons, joins sample metadata, classifies 0N/1N/2N × 3R/4R isoforms, computes per-isoform read counts and `read_count_over_total_sample_read_count`, writes `all_transcripts_processed.csv`. |
| `plot_individual_isoforms_freq_with_stats.py` | Builds the per-tissue ancestry × haplotype heatmaps and box/strip plots, including the `*_all_isoforms_heatmap.png` figures used in the paper. Runs t-test / ANOVA + Tukey HSD with the reference-group structure baked in. |
| `sample_variability_plots.py`                 | Builds the stacked bar plots of per-sample isoform proportions grouped by ancestry-haplotype, one per tissue. |
| `requirements.txt`                            | Python dependencies.                                                                                   |

Each script has a small **CONFIG** block at the top — edit the paths there to
point at your data; there is no separate config file to manage.

## Pipeline

```
BED files + classification .txt files + metadata CSV
                    │
                    ▼
   final_process_all_isoseq_data.py
                    │  writes all_transcripts_processed.csv
                    ▼
   ┌─────────────────────────────────┐
   │ plot_individual_isoforms_freq_  │  → tissue_plots/.../<tissue>/
   │   with_stats.py                 │      ancestry_haplotype_all_isoforms_heatmap.png
   └─────────────────────────────────┘
   ┌─────────────────────────────────┐
   │ sample_variability_plots.py     │  → tissue_plots/sample_variability_plots/
   └─────────────────────────────────┘      sample_variability_<tissue>.png
```

## How to run

1. Edit the **CONFIG** block at the top of each script to point at your input
   directories, metadata file, and desired output location.
2. Install the dependencies:
   ```
   pip install -r requirements.txt
   ```
3. Run the processing step:
   ```
   python final_process_all_isoseq_data.py
   ```
4. Run the plotting scripts:
   ```
   python plot_individual_isoforms_freq_with_stats.py
   python sample_variability_plots.py
   ```

## Inputs

### BED files
PacBio iso-seq BED12 files, one per sample. The processing script extracts
exon coordinates on chr17 (and the H1/H2 alternate contigs:
`chr17_KI270908v1_alt`, `chr17_GL000258v2_alt`) and matches them against the
MAPT exon coordinate tables hard-coded in `final_process_all_isoseq_data.py`
(`h1_exon_mapping`, `h2_exon_mapping`).

### Classification / annotation files
One `*.classification.final.txt` (or `*.clustered.classification.final.txt`)
per sample, in any of the directories listed in `ANNOTATION_DIRS`.
These supply per-isoform Full-Length (`FL` / `Supporting_Reads`) counts and
isoform length. The script auto-detects the on-disk schemas (header at top,
header at bottom, older-format schema).

### Metadata CSV (`METADATA_FILE`)
Required columns:

| Column                  | Description                                                          |
|-------------------------|----------------------------------------------------------------------|
| `File_Name`             | Sample identifier matching the BED filename prefix.                  |
| `HARMONIZED_NAMES`      | Canonical sample name.                                               |
| `Ancestry`              | `EUR`, `AA`, …                                                       |
| `haplotype`             | `H1H1`, `H1H2`, `H1H1_A`, `H2H2`, …                                  |
| `Mutation`              | `WT`. The processing script keeps only `WT` rows.                   |
| `Cell_Type`             | `Astrocytes`, `Forebrain_neurons`, `Brain`, …                        |
| `Brain region_Caudate`  | `x` if caudate, else blank (used to derive `Brain_Caudate` cell type). |
| `Brain region_MFG`      | `x` if MFG, else blank (used to derive `Brain_MFG` cell type).       |

### Reference exon tables
- `H1_transcirpt_name_with_exons.csv`, `H2_transcirpt_name_with_exons.csv`:
  columns `exons`, `source`, `Kat.isoform.definition`, `parsed_exons`,
  `new_transcript_name`. Used to map exon structures to canonical isoform
  names (`0N3R`, `1N3R`, `2N3R_4aL_6`, …).
- `unique_filtered_h1h1_haplotypes_isoforms.csv`,
  `unique_filtered_h2_haplotypes_isoforms.csv`: columns `haplotype`,
  `isoform`, `Transcript_Name`. Used as the haplotype-to-isoform reference.

## Outputs

Everything lands under `OUTPUT_DIR`:

- `all_transcripts_processed.csv` — master per-transcript table (this is the
  input to both plot scripts).
- `isoform_summary.csv`, `harmonized_transcript_exon_inconsistencies.csv` — QC
  summaries.
- `tissue_plots/` — output of the two plotting scripts (heatmaps, box/strip
  plots, stacked bars, Tukey/ANOVA CSVs).

## Notes on the older-format data branch

A subset of the cohort came from an earlier PacBio call set with a different
file schema (header at the bottom; PB-IDs that need a PB→G cross-reference).
The processing script detects these files by checking whether the file path
contains `OLDER_FORMAT_BED_DIR` / `OLDER_FORMAT_ANNOTATION_DIR` (set in the
CONFIG block). If you have no older-format data, leave these set to any
non-matching sentinel string.
