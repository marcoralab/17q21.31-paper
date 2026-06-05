# Figure: GSEA pathway heatmap

Produces `mfg_heatmap_2_added_samples.png` — a heatmap of the significant KEGG
pathways (and their genes) across the haplotype/ancestry comparisons. The MFG
version is the paper figure; the script also builds astro / caudate / neurons
and combined heatmaps.

## How GSEA was computed (`gsea_computation/`)

The pathway tables the heatmap reads are produced by this chain (run per
comparison / per cell type or region — the same scripts were used for each,
only the input paths change):

```
all_genes_result*_with_name_and_type.csv   (DE results from de_analysis/)
  │  1_rank_and_run_gsea.py   ranks genes by the DESeq2 stat and runs
  │                           gseapy.prerank against the MSigDB .gmt sets
  │                           (KEGG / GO-BP / Reactome)
  ▼
gseapy.gene_set.prerank.report.csv   (per comparison: Term, NES, FDR, Lead_genes)
  │  2_pathway_genes_with_fdr.py   attaches per-gene stats to each pathway
  ▼
<label>_pathway_genes_redone_gsea_fdr_q_value_fixed_no_na_padj.csv
  │  3_comprehensive_results.py   aggregates pathways across comparisons,
  │                               classifies effect type  →
  ▼
gsea_comprehensive_results_<tissue>_no_na_padj_flipped_last2.csv
  │  4_filter_significant_genes.py   keeps genes with p < 0.05, counts shared genes
  ▼
<tissue>_filtered_pathways_sig_genes_from_comprehensive_gsea[...].csv
  │  gsea_pathway_heatmap.py
  ▼
mfg_heatmap_2_added_samples.png
```

- **MSigDB gene sets**: `1_rank_and_run_gsea.py` points at `.gmt` files
  (`c2.cp.kegg_legacy`, `c5.go.bp`, `c2.cp.reactome`, v2023.2 human). Download
  these from MSigDB and set the `path/to/msigdb/...` paths. Requires `gseapy`.
- The intermediate file/dir names (`brain/no_na_padj/GSEA_Analysis_...`) and the
  `2_added_samples` suffixes are kept so they match the data layout; edit the
  `path/to/BulkRNASeq/...` prefixes and the relative GSEA-output paths for your
  setup.

## Plotting (`gsea_pathway_heatmap.py`)
Reads the `gsea_comprehensive_results_*` and `*_filtered_pathways_sig_genes_*`
CSVs **by filename from the current working directory** and writes the heatmap
PNGs there. Run it from the directory that holds those files.

```
# per comparison / tissue:
python gsea_computation/1_rank_and_run_gsea.py      # needs: gseapy, pandas
python gsea_computation/2_pathway_genes_with_fdr.py
python gsea_computation/3_comprehensive_results.py
python gsea_computation/4_filter_significant_genes.py
# then:
python gsea_pathway_heatmap.py                      # needs: pandas, numpy, seaborn, matplotlib
```
