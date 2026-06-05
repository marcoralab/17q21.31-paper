"""
Extract structural-variant records from the WhatsHap-phased pbsv VCFs and write
them to a single CSV (Sample, Chromosome, Start, End, SVTYPE).

The CSV is the input to plot_pbsv_variants_ordered_with_clones.py.

Usage:
    python extract_pbsv_sv_to_csv.py \
        --vcf-dir   results/whatshap_output_pbsv \
        --metadata  sample_metadata_with_clones.csv \
        --output    extracted_sv_data_pbsv.csv
"""

import argparse
import os

import pandas as pd
import pysam


def extract_sv_from_vcf(vcf_file: str) -> pd.DataFrame:
    """Return one row per SV record in a single VCF."""
    vcf = pysam.VariantFile(vcf_file)
    records = []
    for record in vcf.fetch():
        if "SVTYPE" not in record.info:
            continue
        records.append(
            {
                "Sample": os.path.basename(vcf_file),
                "Chromosome": record.chrom,
                "Start": record.pos,
                "End": record.stop,
                "SVTYPE": record.info["SVTYPE"],
            }
        )
    return pd.DataFrame(records)


def collect_all_svs(vcf_dir: str) -> pd.DataFrame:
    vcf_files = [
        os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) if f.endswith(".vcf.gz")
    ]
    if not vcf_files:
        raise FileNotFoundError(
            f"No .vcf.gz files found in {vcf_dir}. "
            "Run bgzip_and_index_vcfs.py on the WhatsHap output first."
        )
    return pd.concat(
        [extract_sv_from_vcf(v) for v in vcf_files], ignore_index=True
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf-dir", required=True,
                        help="Directory of phased pbsv .vcf.gz files.")
    parser.add_argument("--metadata", required=True,
                        help="Sample metadata CSV with a 'Sample' column and "
                             "a 'HaplotypePopulation' column.")
    parser.add_argument("--output", required=True,
                        help="Output CSV path.")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    metadata["HaplotypePopulation"] = metadata["HaplotypePopulation"].str.strip()

    sv_data = collect_all_svs(args.vcf_dir)
    sv_data["Sample"] = (
        sv_data["Sample"]
        .str.replace(".all.var.vcf.gz", "", regex=False)
        .str.replace("phased.", "", regex=False)
    )

    combined = pd.merge(metadata, sv_data, on="Sample", how="inner")
    combined["SVTYPE"] = combined["SVTYPE"].fillna("NO_SV")

    sv_data.to_csv(args.output, index=False)
    print(f"Wrote {len(sv_data):,} SV records ({sv_data['Sample'].nunique()} "
          f"samples) to {args.output}")
    print(f"Merged metadata-matched rows: {len(combined):,}")


if __name__ == "__main__":
    main()
