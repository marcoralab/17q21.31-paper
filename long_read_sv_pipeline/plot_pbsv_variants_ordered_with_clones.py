"""
Plot pbsv structural variants across all samples for a given chr17 region,
ordering samples by haplotype population and grouping clones together.

This script produces the figure used in the publication:
    pbsv_variants_region_full_with_clones_ordered_chr17_45307498_46903773.png

Usage:
    python plot_pbsv_variants_ordered_with_clones.py \
        --sv-csv     extracted_sv_data_pbsv.csv \
        --metadata   sample_metadata_with_clones.csv \
        --output-dir output_plots/output_with_clones_ordered
"""

import argparse
import os

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import pandas as pd


# Color scheme used for the publication figure
VARIANT_COLORS = {
    "OTHER": "yellow",
    "SNP":   "pink",
    "INS":   "blue",
    "DEL":   "red",
    "DUP":   "green",
    "BND":   "orange",
    "cnv":   "cyan",
    "INV":   "purple",
}

# Sample-ordering priority by haplotype population
HAPLOTYPE_ORDER = [
    "H1H1_PacBio_EUR",
    "H1H1_PacBio_AFR",
    "H1H1/H1H2_PacBio_AFR",
    "H1H2_PacBio_EUR",
    "H1H2_PacBio_AFR",
    "H2H2_PacBio_EUR",
]

# The 17q21.31 region shown in the publication
REGIONS_OF_INTEREST = [
    {
        "chromosome": "chr17",
        "start": 45307498,
        "end":   46903773,
        "name":  "region_full_with_clones_ordered",
    },
]


def process_region_data(data: pd.DataFrame, region: tuple) -> pd.DataFrame:
    """Subset SVs to the given region and normalize SVTYPE labels."""
    chrom, start, end = region
    region_data = data[
        (data["Chromosome"] == chrom)
        & (data["Start"] >= start)
        & (data["End"] <= end)
    ].copy()

    valid_svtypes = set(VARIANT_COLORS)
    region_data["SVTYPE"] = region_data["SVTYPE"].apply(
        lambda x: x if x in valid_svtypes else "OTHER"
    )
    return region_data


def _clone_sort_key(clone_str: str) -> float:
    """Sort clones numerically (clone1, clone2, ...); non-clones go last."""
    if not clone_str:
        return float("inf")
    try:
        return int(clone_str.replace("clone", ""))
    except ValueError:
        return float("inf")


def ordered_sample_list(metadata: pd.DataFrame) -> list:
    """Build the sample order: by haplotype, clones first within a haplotype."""
    sorted_samples = []
    for haplotype in HAPLOTYPE_ORDER:
        haplotype_samples = metadata[metadata["HaplotypePopulation"] == haplotype]

        clones = haplotype_samples[haplotype_samples["Clone"] != ""].sort_values(
            by="Clone", key=lambda x: x.map(_clone_sort_key)
        )
        sorted_samples.extend(clones["Sample"].tolist())

        non_clones = haplotype_samples[haplotype_samples["Clone"] == ""]
        sorted_samples.extend(non_clones["Sample"].tolist())

    # Reverse so that H1H1_PacBio_EUR is at the top of the plot.
    return sorted_samples[::-1]


def plot_variants(
    data: pd.DataFrame,
    region: tuple,
    output_file: str,
    metadata: pd.DataFrame,
) -> None:
    sorted_samples = ordered_sample_list(metadata)
    sample_to_y = {s: i for i, s in enumerate(sorted_samples)}

    fig, ax = plt.subplots(figsize=(12, 15))

    for sample, y_pos in sample_to_y.items():
        sample_variants = data[data["Sample"] == sample]
        for _, row in sample_variants.iterrows():
            ax.plot(
                row["Start"], y_pos,
                color=VARIANT_COLORS.get(row["SVTYPE"], "gray"),
                marker="o", linestyle="None", markersize=5,
            )

    sample_labels = []
    for sample in sorted_samples:
        clone = metadata.loc[metadata["Sample"] == sample, "Clone"].iloc[0]
        sample_labels.append(f"{sample} ({clone})" if clone else sample)
    plt.yticks(range(len(sorted_samples)), sample_labels, fontsize=6)

    chrom, start, end = region
    ax.set_title(f"Genomic Variants - {chrom}", fontsize=12)
    ax.set_xlabel("Genomic Coordinate", fontsize=10)
    ax.set_ylabel("Samples", fontsize=10)
    ax.legend(
        handles=[
            mlines.Line2D(
                [], [], color=color, marker="o", linestyle="None",
                markersize=5, label=var_type,
            )
            for var_type, color in VARIANT_COLORS.items()
        ],
        loc="upper left", bbox_to_anchor=(1, 1), fontsize=8,
    )
    ax.set_xlim(start, end)
    ax.set_ylim(-1, len(sorted_samples))

    plt.subplots_adjust(right=0.75)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-csv", required=True,
                        help="CSV from extract_pbsv_sv_to_csv.py.")
    parser.add_argument("--metadata", required=True,
                        help="Sample metadata CSV with Sample, "
                             "HaplotypePopulation, and Clone columns.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory to write the figure(s) into.")
    args = parser.parse_args()

    print("Loading data ...")
    metadata = pd.read_csv(args.metadata)
    pbsv_data = pd.read_csv(args.sv_csv)

    metadata["Sample"] = metadata["Sample"].str.strip("-")
    metadata["HaplotypePopulation"] = metadata["HaplotypePopulation"].str.strip()
    metadata["Clone"] = metadata["Clone"].fillna("")

    pbsv_data["Sample"] = (
        pbsv_data["Sample"]
        .str.replace(".all.var.vcf.gz", "", regex=False)
        .str.replace("phased.", "", regex=False)
    )

    os.makedirs(args.output_dir, exist_ok=True)

    for region in REGIONS_OF_INTEREST:
        chrom = region["chromosome"]
        start = region["start"]
        end = region["end"]
        name = region["name"]

        region_data = process_region_data(pbsv_data, (chrom, start, end))
        output_file = os.path.join(
            args.output_dir, f"pbsv_variants_{name}_{chrom}_{start}_{end}.png"
        )
        plot_variants(region_data, (chrom, start, end), output_file, metadata)
        print(f"Saved {chrom}:{start}-{end} -> {output_file}")


if __name__ == "__main__":
    main()
