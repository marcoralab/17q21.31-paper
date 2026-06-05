# Bulk RNA-seq: 17q21 inversion / MAPT haplotype × ancestry analysis

Code for the bulk RNA-seq analyses of *MAPT* H1/H2 haplotype and genetic
ancestry, including the differential-expression pipeline and the figures built
from it. The same pipeline was applied across the iPSC-derived cell types
(astrocytes, neurons, microglia) and the post-mortem brain regions (DLPFC, MFG,
caudate); only the input/output paths differ per dataset.

The repository is organized **by figure**. Two folders are shared upstream
stages (`alignment/`, `de_analysis/`); the rest each reproduce one published
figure.

```
bulkrnaseq_mapt_inversion/
├── alignment/                STAR → Picard → RSEM (FASTQ → gene counts)
├── de_analysis/              DESeq2 + CQN differential expression, gene annotation
├── figure_inversion_genes/   shared-genes dot plot + volcano plots (inversion locus)
├── figure_gsea_heatmap/      GSEA pathway heatmap (+ gsea_computation/ to run GSEA)
└── figure_venn_shared_genes/ Venn diagram of shared DE genes (any cell type/region)
```

> **Data is not bundled.** Every script reads its inputs from paths set in a
> CONFIG block (or clearly-marked `path/to/...` placeholders) at the top of the
> file. Point them at your downloaded data. The starting data is provided
> separately.

## Pipeline overview

```
FASTQ
  │  alignment/  (STAR → Picard QC → RSEM)
  ▼
per-sample *.genes.results
  │  de_analysis/differential_expression_deseq2.R   (DESeq2 + CQN, one run per contrast)
  ▼
all_genes_result*.csv   (gene_id, log2FoldChange, pvalue, padj, ...)
  │  de_analysis/make_gene_id_mapping.py + annotate_de_results.py
  ▼
*_with_name_and_type.csv   (adds Gene_Symbol, Gene_Type)  ← the figure inputs
  │
  ├── figure_inversion_genes/   filter_to_inversion_genes.R → shared_genes_dotplot.py + volcano_plots.py
  ├── figure_gsea_heatmap/      gsea_computation/ (rank → GSEA → pathway tables) → gsea_pathway_heatmap.py
  └── figure_venn_shared_genes/ venn_shared_genes.R
```

## Shared stages

### `alignment/`
STAR alignment → Picard RNA-seq QC metrics → RSEM quantification. Each stage has
a `*_create_params.sh` that scans a FASTQ/BAM directory and writes a parameter
file, and a job script that consumes it. The same scripts were run per cell type
(astrocytes / neurons / microglia) and per brain region (DLPFC / MFG / caudate);
only the input/output directories differ. Edit the path variables at the top of
each script. Requires `STAR`, `picard`, `RSEM`, `samtools`.

### `de_analysis/`
- `differential_expression_deseq2.R` — reads the per-sample RSEM `.genes.results`,
  filters/normalizes (CQN with Ensembl-111 GC content), and runs DESeq2 for the
  haplotype / ancestry contrasts, writing `all_genes_result*.csv` per contrast.
- `make_gene_id_mapping.py` — parses the GENCODE GTF into a
  `gene_id → (symbol, type)` table.
- `annotate_de_results.py` — uses that table to add `Gene_Symbol`/`Gene_Type`
  to a DESeq2 result CSV, producing the `*_with_name_and_type.csv` the figure
  scripts read.

## Figures

| Folder | Figure | Script(s) |
|--------|--------|-----------|
| `figure_inversion_genes/` | `shared_genes_downregulated_MFG_..._2_MFG_samples.png` and the inversion-locus volcano plots | `filter_to_inversion_genes.R` → `shared_genes_dotplot.py`, `volcano_plots.py` |
| `figure_gsea_heatmap/` | `mfg_heatmap_2_added_samples.png` | `gsea_computation/` (1→4) → `gsea_pathway_heatmap.py` |
| `figure_venn_shared_genes/` | `Downregulated_Genes_Venn_<label>.png` (e.g. `_neurons`) | `venn_shared_genes.R` |

Each figure folder has its own README with the exact inputs and run command.

## Notes

- The differential-expression script (`differential_expression_deseq2.R`) was
  cleaned for paths and de-cruft (removed REPL/conda lines and the cluster proxy
  settings) but its internal variable names were **left as-is** to guarantee the
  results are byte-identical to the published analysis. It is run once per
  contrast (the contrast is selected inside the script).
- The `*_2_added_samples` / `2_added_MFG_samples` naming throughout refers to two
  MFG samples added to that cohort late; it is kept so filenames match the data.
