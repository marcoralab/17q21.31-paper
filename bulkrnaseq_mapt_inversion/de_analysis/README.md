# Differential expression & gene annotation

Turns per-sample RSEM gene counts into the annotated differential-expression
tables the figure scripts consume.

| Script | Role |
|--------|------|
| `differential_expression_deseq2.R` | Reads the RSEM `.genes.results`, filters to chr1–Y genes, normalizes with CQN (Ensembl-111 GC content), and runs DESeq2 for the haplotype / ancestry contrasts. Writes `all_genes_result*.csv` (and `.rda` checkpoints). Edit the CONFIG paths and `setwd()` target at the top; run once per contrast. |
| `make_gene_id_mapping.py` | Parses the GENCODE GTF → `gene_id_to_symbol_and_type_mapping.csv` (`ENSEMBL_ID, Gene_Symbol, Gene_Type`). |
| `annotate_de_results.py` | Adds `Gene_Symbol` / `Gene_Type` to a DESeq2 result CSV using that mapping → `*_with_name_and_type.csv` (the form the figure scripts read). List the CSVs to annotate in the CONFIG `input_files`. |

## Run
```
Rscript differential_expression_deseq2.R   # needs: DESeq2, cqn, biomaRt, limma, edgeR, variancePartition, corrplot
python make_gene_id_mapping.py
python annotate_de_results.py
```

> Note: `differential_expression_deseq2.R` had its paths replaced with
> placeholders and its REPL/conda/proxy lines removed, but the internal variable
> names were intentionally left unchanged so results stay identical to the
> published analysis.
