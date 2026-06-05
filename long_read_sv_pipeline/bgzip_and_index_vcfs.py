"""
bgzip + tabix-index every .vcf file in a directory.

The downstream extraction script (extract_pbsv_sv_to_csv.py) expects bgzipped
.vcf.gz files, so run this once on the WhatsHap-phased pbsv output directory
after the Snakemake workflow finishes.

Usage:
    python bgzip_and_index_vcfs.py <vcf_directory>
"""

import argparse
import os
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "vcf_dir",
        help="Directory containing .vcf files to bgzip and tabix-index.",
    )
    args = parser.parse_args()

    vcf_files = [f for f in os.listdir(args.vcf_dir) if f.endswith(".vcf")]
    for vcf_file in vcf_files:
        vcf_path = os.path.join(args.vcf_dir, vcf_file)
        vcf_gz_path = vcf_path + ".gz"

        print(f"Compressing {vcf_file} ...")
        with open(vcf_gz_path, "wb") as fh:
            subprocess.run(["bgzip", "-c", vcf_path], stdout=fh, check=True)

        print(f"Indexing {vcf_file}.gz ...")
        subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

    print("Compression and indexing completed.")


if __name__ == "__main__":
    main()
