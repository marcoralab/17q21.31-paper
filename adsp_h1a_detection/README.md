# ADSP MAPT H1/H2 haplotype × Alzheimer's disease risk

Code to test the association between *MAPT* H1/H2 haplotype and Alzheimer's
disease (AD) risk in the ADSP r4 WGS cohort, stratified by genetic ancestry,
using a GENESIS mixed-model two-step (Cholesky-residual) analysis.

The pipeline reproduces the publication figure:

```
results/plots/haplotype_forest_plot_cholesky_filtered.png
```

a forest plot of haplotype coefficient estimates (95% CI, reference = H1H1)
for the African, European, and overall samples.

## Pipeline

```
ADSP chr17 VCF
      │  genotype_adsp_data.py   (bcftools extract of 2 MAPT tag SNPs)
      ▼
results/mapt_genotypes_v10.csv
      │  agree_disagree_evaluation.py   (tag-SNP agreement + haplotype call)
      ▼
results/mapt_genotypes_v10_with_agreement.csv
      │  merge_with_phenotype_file.py   (+ phenotype TSV, ancestry TSV, covariate/PC TSV)
      ▼
results/mapt_genotypes_with_covar_merged.csv
      │  model_genesis_haplotype_ad.R   (GENESIS null model → Cholesky 2-step → forest plot)
      ▼
results/plots/haplotype_forest_plot_cholesky_filtered.png
```

Each script writes into `results/` relative to the working directory, so run
them from this folder, in order.

## Files

| File                            | Purpose                                                                                             |
|---------------------------------|-----------------------------------------------------------------------------------------------------|
| `genotype_adsp_data.py`         | Extracts the two MAPT H1/H2 tag SNPs (`rs8070723`, `rs1052553`) from the ADSP chr17 VCF via `bcftools`. |
| `agree_disagree_evaluation.py`  | Flags tag-SNP agreement and derives a per-sample `Haplotype` / `sample_name_cut`.                   |
| `merge_with_phenotype_file.py`  | Merges genotypes with the phenotype, global-ancestry, and PCA-covariate tables into the modeling input. |
| `model_genesis_haplotype_ad.R`  | Fits the GENESIS GLMM null model (kinship + sequencing random effects), runs the Cholesky-residual second step per ancestry, and draws the forest plot. |
| `requirements.txt`              | Python dependencies (plus external `bcftools`).                                                     |

Each script has a **CONFIG** block at the top — edit the paths there to point
at your data.

## How to run

```
# 1. Extract MAPT tag-SNP genotypes (edit VCF_FILE first)
python genotype_adsp_data.py

# 2. Call haplotype agreement
python agree_disagree_evaluation.py

# 3. Merge phenotype + ancestry + covariates (edit the three *_TSV paths first)
python merge_with_phenotype_file.py

# 4. Fit the model and draw the forest plot (edit grm_path first)
Rscript model_genesis_haplotype_ad.R
```

## Inputs

All inputs except the ADSP genotypes come from controlled-access / pipeline
files; point the CONFIG paths at your equivalents.

- **`VCF_FILE`** — bgzipped, tabix-indexed chr17 VCF from the ADSP r4 WGS
  release (NIAGADS/dbGaP controlled access).
- **`PHENOTYPE_TSV`** — sample annotation keyed on `sample.id`, supplying at
  least `diagnosis`, `sex`, `apoe`, `sequencing_center`, `sequencing_platform`.
- **`ANCESTRY_TSV`** — global-ancestry / admixture estimates keyed on `IID`,
  including the ancestry fractions and `admixture_super_pop_max`.
- **`COVARIATE_TSV`** — sample annotation with the PCA covariate columns
  `V1..Vn`, keyed on `sample.id`.
- **`grm_path`** — a sample × sample genetic relatedness matrix saved as `.rds`
  (e.g. a PC-Relate / GCTA kinship matrix) with sample IDs as row/column names.

### Columns expected by the R script in the merged CSV
`sample.id`, `diagnosis`, `sex`, `apoe`, `admixture_super_pop_max`,
`Agreement`, `Haplotype`, `rs8070723`, `sequencing_center`,
`sequencing_platform`, and the PC columns `V1..V5`.

## Analysis notes

- **Haplotype coding.** Samples are grouped by tag-SNP genotype into H1H1 /
  H1H2 / H2H2; samples whose two tag SNPs disagree are treated as an
  `H1_A Carrier` group. H1H1 is the reference level.
- **Two-step model.** Step 1 fits a GENESIS GLMM for AD with the kinship matrix
  and sequencing platform as random effects and base covariates + PCs as fixed
  effects, *without* haplotype. Step 2 regresses the Cholesky residuals on
  haplotype, giving the per-haplotype coefficient estimates plotted.
- The `_filtered` forest plot drops the African H2H2 estimate (too few samples
  for a stable estimate); the unfiltered `_cholesky_only` plot is also written.

## Dependencies

- Python ≥ 3.8 with `pandas` (see `requirements.txt`), plus `bcftools` on PATH.
- R with: `GENESIS`, `Biobase`, `Matrix`, `lme4`, `tidyverse` (`dplyr`,
  `ggplot2`, `tidyr`), `magrittr`.
