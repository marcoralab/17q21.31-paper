# Figure: shared DE-gene Venn diagram

A 5-way Venn of the significantly **downregulated** genes shared across five
haplotype/ancestry comparisons. The script also writes the upregulated and
combined Venns and the set-intersection tables.

This was used for several datasets (cell types and brain regions). Set `label`
in the CONFIG block to the dataset you are plotting; `label = "neurons"`
reproduces the paper figure `Downregulated_Genes_Venn_neurons.png`.

## Script
`venn_shared_genes.R`. In the CONFIG block set `label`, the five input CSV
paths, and `output_dir`.

## Inputs
Five DESeq2 result tables (`*_with_name_and_type.csv` from `de_analysis/`) for
the chosen dataset, one per comparison:
`H1H1A_AA vs H1H1_AA`, `H1H2_AA vs H1H2_EUR`, `H1H2_EUR vs H1H1_EUR`,
`H1H2_AA vs H1H1_AA`, `H1H1_AA vs H1H1_EUR`.
Each must contain `Gene_Symbol`, `log2FoldChange`, `pvalue`, `padj`.
Significance threshold: `pvalue < 0.05` and `|log2FoldChange| > 1`.

## Run
```
Rscript venn_shared_genes.R     # needs: VennDiagram, dplyr, grid
```
Outputs `Downregulated_Genes_Venn_<label>.png` (+ up / combined Venns and
intersection CSVs) in `output_dir`.

> Cleaned from the original, which loaded several datasets in sequence and was
> run interactively. This version keeps the identical filtering logic but is
> parameterized by `label` so each dataset's figure is produced directly. The
> unrelated overlap-heatmap tail of the original (extra packages + syntax
> errors) was dropped.
