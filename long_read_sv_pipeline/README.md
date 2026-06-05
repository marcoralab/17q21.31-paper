# Long-read SV calling, phasing, and 17q21.31 plotting

This directory contains the long-read code used to produce the publication
figure
`pbsv_variants_region_full_with_clones_ordered_chr17_45307498_46903773.png`
(chr17:45,307,498–46,903,773, samples ordered by haplotype population with
clones grouped).

## Workflow

1. **Align + call + phase (Snakemake).**
   `Snakefile` aligns each PacBio HiFi FASTQ to GRCh38 with `pbmm2`, calls
   structural variants with `pbsv discover` / `pbsv call`, and phases the
   resulting VCFs with `WhatsHap`.

   ```
   snakemake --use-conda --cores <N>
   ```

   Conda environments are under `envs/`.

2. **Compress and index the phased pbsv VCFs.**

   ```
   python bgzip_and_index_vcfs.py results/whatshap_output_pbsv
   ```

3. **Extract SV records to a single CSV.**

   ```
   python extract_pbsv_sv_to_csv.py \
       --vcf-dir   results/whatshap_output_pbsv \
       --metadata  sample_metadata_with_clones.csv \
       --output    extracted_sv_data_pbsv.csv
   ```

4. **Produce the publication figure.**

   ```
   python plot_pbsv_variants_ordered_with_clones.py \
       --sv-csv     extracted_sv_data_pbsv.csv \
       --metadata   sample_metadata_with_clones.csv \
       --output-dir output_plots/output_with_clones_ordered
   ```

   This writes
   `output_plots/output_with_clones_ordered/pbsv_variants_region_full_with_clones_ordered_chr17_45307498_46903773.png`.

## Sample metadata (`sample_metadata_with_clones.csv`)

Required columns:

| Column                | Description                                                  |
|-----------------------|--------------------------------------------------------------|
| `Sample`              | Sample identifier; must match the FASTQ keys in `config.yaml`. |
| `HaplotypePopulation` | e.g. `H1H1_PacBio_EUR`, `H1H2_PacBio_AFR`, …                 |
| `Haplotype`           | Inferred 17q21.31 haplotype (`H1H1_PacBio`, …).              |
| `Population`          | `EUR`, `AFR`, …                                              |
| `rs8070723`           | Allele at rs8070723.                                         |
| `rs1052553`           | Allele at rs1052553.                                         |
| `Clone`               | Clone label (`clone1`, `clone2`, …) or empty for non-clones. |

## Inputs and outputs

- **Reference**: GRCh38 primary assembly (`GRCh38.primary_assembly.genome.fa`
  and its `.mmi` index). Paths set in `config.yaml`.
- **Input**: one PacBio HiFi FASTQ per sample, listed under `samples:` in
  `config.yaml`.
- **Key output**: phased pbsv VCFs in
  `results/whatshap_output_pbsv/phased.<sample>.all.var.vcf` and the figure
  PNG above.

## Software versions

See `envs/sv.yaml` and `envs/whatshap.yaml`. Pinned versions used for the
publication: pbmm2 1.12.0, pbsv 2.2.0, WhatsHap 2.0.
