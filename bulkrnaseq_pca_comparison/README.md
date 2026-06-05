# Bulk RNA-seq PCA: iPSC-derived cells vs. public datasets

PCA comparing this study's CQN-normalized iPSC-derived astrocyte and neuron
bulk RNA-seq against five public reference datasets, to show where the
iPSC-derived cells sit relative to published iPSC and primary cell data.

Reproduces the publication figure:

```
results/plots/PCs_plot_five_shapes_no_circle.png
```

a PC1 vs PC2 scatter where **color = sample group** (haplotype × ancestry ×
cell type, plus each external dataset) and **shape = dataset of origin**
(five non-circle shapes: this study's iPSC astrocytes & neurons, Bowles 2024,
Leng 2022, TCW 2022). The script also writes a sister plot
`PCs_plot_with_ancestry_lines_with_leng.png` (color only) from the same PCA.

## File

| File | Purpose |
|------|---------|
| `pca_combined_rnaseq_comparison.R` | Loads this study's normalized expression + the five public datasets, takes the shared genes, runs PCA on the combined log-expression matrix, and draws the PCA plots. |

Edit the **CONFIG** block at the top of the script to point at your files.

## How to run

```
Rscript pca_combined_rnaseq_comparison.R
```

Outputs are written to `results/plots/`.

## Inputs

### This study's normalized expression (`*_RDA`)
Three `.rda` files (one per cell type) produced by the gene-expression
filtering / CQN-normalization step of the bulk RNA-seq pipeline. Each restores:

- `datExpr_astro` / `datExpr_neuron` / `datExpr` — genes (Ensembl IDs, rows) ×
  samples (columns), CQN-normalized log-expression.
- `datMeta_astro` / `datMeta_neuron` / `datMeta` — sample metadata with at
  least `File.Name_1`, `haplotype` (`H1H1`, `H1H2`, `H1H1_A`), and `Ancestry`
  (`EUR`, `AA`).

> The microglia file is loaded but its samples are **not** plotted — it is kept
> only so the cross-dataset common-gene set (and therefore the PCA) matches the
> published figure.

### Public comparison datasets (download from GEO)

| CONFIG variable  | Source | File |
|------------------|--------|------|
| `GSE260517_FILE` | Bowles 2024 (GSE260517) | STAR gene counts (`GSE260517_star_gene.txt`) |
| `GSE190185_FILE` | TCW 2022 (GSE190185)    | featureCounts (`GSE190185_isogenicAPOE_featureCounts.txt`) |
| `GSE73721_FILE`  | Zhang 2016 (GSE73721)   | gene table with Ensembl IDs (`.csv`) |
| `GSE182307_FILE` | Leng 2022 (GSE182307)   | gene-level counts (`GSE182307_geneLevelCounts_allSamples.csv`) |
| `GSE97904_FILE`  | GSE97904                | expression matrix (`GSE97904_expression.tsv`) |

The specific sample columns used from each dataset are selected explicitly in
the script (see the `selected_samples_*` vectors).

## Method

1. Ensembl gene IDs are stripped of version suffixes across all datasets.
2. Each external counts matrix is normalized with edgeR (`calcNormFactors` →
   `cpm(log = TRUE)`); this study's matrices are already CQN-normalized.
3. The datasets are intersected to their common genes and column-bound into one
   matrix, which is gene-centered (`scale(scale = FALSE)`) and run through
   `prcomp`.
4. Samples are filtered to the plotted groups, the PCA is recomputed, and PC1
   vs PC2 is plotted with per-group colors and per-dataset shapes.

## Dependencies

R with: `tidyverse` (`dplyr`, `ggplot2`, `tibble`), `limma`, `edgeR`, `sva`,
`biomaRt`, `RColorBrewer`.
