# Figure: inversion-locus shared-genes dot plot + volcano plots

Two figures built from the 17q21-inversion-locus genes in the MFG comparisons.

## Scripts (run in order)

1. **`filter_to_inversion_genes.R`** — subsets each differential-expression
   result to genes in the inversion locus (using the NCBI inversion gene list
   and a GENCODE GTF for id→name mapping), re-computes BH-adjusted p-values, and
   writes `inversion_<comparison>.csv`. Set `de_dir` to the folder of annotated
   DE results (named `<comparison>.csv`, using the names in the `comparisons`
   list), and the GTF / gene-list paths, in the CONFIG block.

2. **`shared_genes_dotplot.py`** — reads the inversion-filtered MFG comparison
   CSVs and draws the binary presence/absence dot plot of significant genes
   across comparisons →
   `shared_genes_downregulated_MFG_with_flipped_geneslast2_added_2_MFG_samples.png`
   (and the upregulated variant).

3. **`volcano_plots.py`** — reads the same inversion-filtered MFG CSVs and draws
   one volcano plot per comparison (pvalue and padj variants).

For `shared_genes_dotplot.py` and `volcano_plots.py`, set `base_dir` / `data_dir`
to the folder of inversion-filtered CSVs and `output_dir` to where plots go.

## Inputs
The MFG `inversion_*_brain_mfg*.csv` files produced by
`filter_to_inversion_genes.R`. Each needs columns `gene_name`,
`log2FoldChange`, `pvalue`, `padj`.

## Run
```
Rscript filter_to_inversion_genes.R
python shared_genes_dotplot.py
python volcano_plots.py     # needs: pandas, numpy, matplotlib, adjustText
```
