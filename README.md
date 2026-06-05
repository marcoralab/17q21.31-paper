# MAPT / 17q21.31 inversion — analysis code

This repository includes the code for long-read sequencing, iso-seq, bulk RNA-seq, and the ADSP analysis.

## Repository structure

| Folder | Analysis | Key figure(s) |
|--------|----------|---------------|
| [`long_read_sv_pipeline/`](long_read_sv_pipeline) | PacBio HiFi long-read alignment, `pbsv` structural-variant calling and `WhatsHap` phasing across the 17q21.31 region (Snakemake). | `pbsv_variants_region_full_with_clones_ordered_*.png` |
| [`isoseq_analysis/`](isoseq_analysis) | PacBio Iso-Seq exon-level *MAPT* isoform calling (0N/1N/2N × 3R/4R) and tissue-stratified isoform-usage plots. | ancestry-haplotype isoform heatmaps; per-sample isoform-proportion bars |
| [`bulkrnaseq_mapt_inversion/`](bulkrnaseq_mapt_inversion) | Bulk RNA-seq (iPSC-derived cell types + post-mortem brain regions): STAR→Picard→RSEM, DESeq2/CQN differential expression, GSEA, and the inversion-locus / pathway / Venn figures. | inversion shared-gene dot plot, volcano plots, GSEA pathway heatmap, shared-DE-gene Venn |
| [`bulkrnaseq_pca_comparison/`](bulkrnaseq_pca_comparison) | PCA comparing the study's iPSC-derived astrocyte/neuron bulk RNA-seq against five public reference datasets. | `PCs_plot_five_shapes_no_circle.png` |
| [`adsp_h1a_detection/`](adsp_h1a_detection) | ADSP r4 WGS: *MAPT* H1/H2 tag-SNP genotyping and a GENESIS mixed-model (Cholesky two-step) test of haplotype × Alzheimer's-disease risk, stratified by ancestry. | `haplotype_forest_plot_cholesky_filtered.png` |

## How to use this code

- **Each folder is independent.** Start with that folder's `README.md` — it lists
  the scripts, the run order, the inputs, and the dependencies for that analysis.
- **Paths are not hard-coded.** Every script has a `CONFIG` block (or clearly
  marked `path/to/...` placeholders) at the top.
- **Data is not included here.** 
- **Languages / tools vary by analysis** — Python, R/Bioconductor, and
  cluster tools (STAR, RSEM, Picard, pbsv, WhatsHap, GENESIS, gseapy). Each
  folder's README lists what it needs.

## Citation

If you use this code, please cite the paper.
