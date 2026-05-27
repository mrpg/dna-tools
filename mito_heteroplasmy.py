#!/usr/bin/env python3
"""
Mitochondrial heteroplasmy landscape analysis.

Plots heteroplasmy frequencies across the mitogenome, identifies age-related
hotspots, and screens against MITOMAP pathogenic mutations.

Usage:
    python mito_heteroplasmy.py <mito.vcf.gz> <output_dir> [--sample-name NAME]

Input VCF must have FORMAT fields: GT, DP, HF (heteroplasmy frequency).
Requires: matplotlib, numpy, bcftools in PATH.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# Mitochondrial genome annotation (rCRS)
REGIONS = [
    (1, 576, "D-loop", "regulatory"),
    (577, 647, "tRNA-Phe", "trna"),
    (648, 1601, "12S rRNA", "rrna"),
    (1602, 1670, "tRNA-Val", "trna"),
    (1671, 3229, "16S rRNA", "rrna"),
    (3230, 3304, "tRNA-Leu", "trna"),
    (3307, 4262, "ND1", "complex_i"),
    (4470, 5511, "ND2", "complex_i"),
    (5904, 7445, "COX1", "complex_iv"),
    (7586, 8269, "COX2", "complex_iv"),
    (8366, 8572, "ATP8", "complex_v"),
    (8527, 9207, "ATP6", "complex_v"),
    (9207, 9990, "COX3", "complex_iv"),
    (10059, 10404, "ND3", "complex_i"),
    (10470, 10766, "ND4L", "complex_i"),
    (10760, 12137, "ND4", "complex_i"),
    (12337, 14148, "ND5", "complex_i"),
    (14149, 14673, "ND6", "complex_i"),
    (14747, 15887, "CytB", "complex_iii"),
    (16024, 16569, "D-loop", "regulatory"),
]

REGION_COLORS = {
    "regulatory": "#E8D44D",
    "trna": "#E8773D",
    "rrna": "#7DB3E0",
    "complex_i": "#6BBF6B",
    "complex_iii": "#C97DBF",
    "complex_iv": "#E05D5D",
    "complex_v": "#5DC1C1",
}

# Known age-related heteroplasmy hotspots (MITOMAP / literature)
AGE_HOTSPOTS = {
    150: "m.150C>T \u2014 age-associated, D-loop",
    152: "m.152T>C \u2014 age-associated, D-loop",
    189: "m.189A>G \u2014 age-associated, D-loop",
    195: "m.195T>C \u2014 age-associated, D-loop",
    204: "m.204T>C \u2014 age-associated, D-loop",
    408: "m.408T>A \u2014 age-associated, D-loop",
    414: "m.414T>G \u2014 age-associated, D-loop",
    16093: "m.16093T>C \u2014 age-associated, D-loop",
}

# MITOMAP confirmed pathogenic mutations (subset of most well-established)
PATHOGENIC = {
    583: ("G>A", "MELAS", "tRNA-Phe"),
    1555: ("A>G", "Aminoglycoside deafness", "12S rRNA"),
    1624: ("C>T", "MELAS", "tRNA-Val"),
    3243: ("A>G", "MELAS / Diabetes", "tRNA-Leu(UUR)"),
    3256: ("C>T", "MELAS", "tRNA-Leu(UUR)"),
    3260: ("A>G", "Cardiomyopathy", "tRNA-Leu(UUR)"),
    3271: ("T>C", "MELAS", "tRNA-Leu(UUR)"),
    3460: ("G>A", "LHON", "ND1"),
    4171: ("C>A", "LHON", "ND1"),
    5537: ("insT", "Encephalomyopathy", "tRNA-Trp"),
    5545: ("C>T", "Encephalopathy", "tRNA-Trp"),
    7445: ("A>G", "Deafness", "COX1/tRNA-Ser"),
    7472: ("insC", "Deafness + Neurological", "tRNA-Ser(UCN)"),
    8344: ("A>G", "MERRF", "tRNA-Lys"),
    8356: ("T>C", "MERRF", "tRNA-Lys"),
    8363: ("G>A", "Cardiomyopathy", "tRNA-Lys"),
    8993: ("T>C", "NARP / Leigh syndrome", "ATP6"),
    9176: ("T>C", "Leigh syndrome", "ATP6"),
    10158: ("T>C", "Leigh syndrome", "ND3"),
    10191: ("T>C", "Leigh syndrome", "ND3"),
    10197: ("G>A", "Leigh / Dystonia", "ND3"),
    11778: ("G>A", "LHON", "ND4"),
    13513: ("G>A", "MELAS / Leigh", "ND5"),
    14459: ("G>A", "LHON + Dystonia", "ND6"),
    14484: ("T>C", "LHON", "ND6"),
    14487: ("T>C", "Leigh / Dystonia", "ND6"),
    14709: ("T>C", "Myopathy / Diabetes", "tRNA-Glu"),
    15990: ("C>T", "Exercise intolerance", "tRNA-Pro"),
}


def matches_pathogenic_change(variant, expected_change):
    """Return True when a VCF allele matches the MITOMAP change exactly."""
    ref = variant["ref"]
    alt = variant["alt"]

    if ">" in expected_change:
        expected_ref, expected_alt = expected_change.split(">", 1)
        return ref == expected_ref and alt == expected_alt

    if expected_change.startswith("ins"):
        inserted = expected_change.removeprefix("ins")
        return alt.startswith(ref) and alt[len(ref) :] == inserted

    return False


def extract_variants_from_vcf(vcf_path):
    """Extract variants from a mito VCF using bcftools."""
    result = subprocess.run(
        ["bcftools", "query", "-f", r"%POS\t%REF\t%ALT\t[%HF]\t[%DP]\n", str(vcf_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    variants = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        pos = int(parts[0])
        ref = parts[1]
        alts = parts[2].split(",")
        hfs = parts[3].split(",")
        dp = int(parts[4].split(",")[0])
        for alt, hf in zip(alts, hfs):
            variants.append(
                {
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "hf": float(hf),
                    "dp": dp,
                }
            )
    return variants


def get_sample_name_from_vcf(vcf_path):
    """Extract sample name from VCF header."""
    result = subprocess.run(
        ["bcftools", "query", "-l", str(vcf_path)],
        capture_output=True,
        text=True,
    )
    names = result.stdout.strip().split("\n")
    return names[0] if names else "Unknown"


def get_region_type(pos):
    for start, end, name, rtype in REGIONS:
        if start <= pos <= end:
            return name, rtype
    return "Intergenic", "regulatory"


def plot_landscape(variants, out_dir):
    """Create the heteroplasmy landscape plot."""
    fig, axes = plt.subplots(
        3, 1, figsize=(16, 12), gridspec_kw={"height_ratios": [1, 4, 2]}
    )

    # --- Top panel: Gene map ---
    ax_genes = axes[0]
    ax_genes.set_xlim(0, 16569)
    ax_genes.set_ylim(0, 1)
    ax_genes.set_yticks([])
    ax_genes.set_title(
        "Mitochondrial Genome Map (16,569 bp)", fontsize=13, fontweight="bold"
    )

    for start, end, name, rtype in REGIONS:
        color = REGION_COLORS.get(rtype, "#CCCCCC")
        ax_genes.barh(
            0.5,
            end - start,
            left=start,
            height=0.6,
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        if (end - start) > 500:
            ax_genes.text(
                (start + end) / 2, 0.5, name, ha="center", va="center", fontsize=6
            )

    legend_patches = [
        mpatches.Patch(color=REGION_COLORS["complex_i"], label="Complex I (ND1-6)"),
        mpatches.Patch(color=REGION_COLORS["complex_iii"], label="Complex III (CytB)"),
        mpatches.Patch(color=REGION_COLORS["complex_iv"], label="Complex IV (COX1-3)"),
        mpatches.Patch(color=REGION_COLORS["complex_v"], label="Complex V (ATP6/8)"),
        mpatches.Patch(color=REGION_COLORS["rrna"], label="rRNA"),
        mpatches.Patch(color=REGION_COLORS["trna"], label="tRNA"),
        mpatches.Patch(color=REGION_COLORS["regulatory"], label="D-loop"),
    ]
    ax_genes.legend(handles=legend_patches, loc="upper right", fontsize=7, ncol=4)

    # --- Middle panel: Heteroplasmy frequencies ---
    ax_hf = axes[1]
    low_level = [v for v in variants if v["hf"] < 0.95]
    homoplasmic = [v for v in variants if v["hf"] >= 0.95]

    colors_low = [
        REGION_COLORS.get(get_region_type(v["pos"])[1], "#999999") for v in low_level
    ]
    ax_hf.scatter(
        [v["pos"] for v in low_level],
        [v["hf"] for v in low_level],
        c=colors_low,
        s=8,
        alpha=0.6,
        edgecolors="none",
    )

    for v in homoplasmic:
        ax_hf.axvline(
            v["pos"], color="#333333", alpha=0.3, linewidth=0.5, linestyle="--"
        )

    # Annotate age-related hotspots
    for v in low_level:
        if v["pos"] in AGE_HOTSPOTS and v["hf"] > 0.005:
            ax_hf.annotate(
                f"m.{v['pos']}\nHF={v['hf']:.3f}",
                (v["pos"], v["hf"]),
                textcoords="offset points",
                xytext=(10, 10),
                fontsize=7,
                color="#CC3333",
                arrowprops=dict(arrowstyle="->", color="#CC3333", lw=0.8),
            )

    # Annotate pathogenic hits
    for v in low_level:
        if v["pos"] in PATHOGENIC:
            ax_hf.annotate(
                f"m.{v['pos']}\n\u26a0 {PATHOGENIC[v['pos']][1]}",
                (v["pos"], v["hf"]),
                textcoords="offset points",
                xytext=(10, -15),
                fontsize=7,
                color="#FF6600",
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#FF6600", lw=1),
            )

    hf_low = [v["hf"] for v in low_level]
    ax_hf.set_xlim(0, 16569)
    ax_hf.set_ylim(0, max(hf_low) * 1.15 if hf_low else 0.05)
    ax_hf.set_ylabel("Heteroplasmy Frequency (HF)", fontsize=11)
    ax_hf.set_title(
        "Heteroplasmy Landscape (sub-clonal variants, HF < 0.95)",
        fontsize=13,
        fontweight="bold",
    )
    ax_hf.axhline(
        0.01, color="#999999", linestyle=":", linewidth=0.8, label="1% threshold"
    )
    ax_hf.legend(fontsize=8)

    # --- Bottom panel: Read depth ---
    ax_dp = axes[2]
    unique_pos = sorted(set(v["pos"] for v in variants))
    pos_to_dp = {v["pos"]: v["dp"] for v in variants}
    ax_dp.fill_between(
        unique_pos, [pos_to_dp[p] for p in unique_pos], alpha=0.4, color="#4477AA"
    )
    ax_dp.scatter(
        [v["pos"] for v in variants],
        [v["dp"] for v in variants],
        s=2,
        alpha=0.3,
        color="#4477AA",
    )
    ax_dp.set_xlim(0, 16569)
    ax_dp.set_ylabel("Read Depth", fontsize=11)
    ax_dp.set_xlabel("Mitochondrial Position (bp)", fontsize=11)
    ax_dp.set_title("Sequencing Depth", fontsize=13, fontweight="bold")

    plt.tight_layout()
    plt.savefig(out_dir / "heteroplasmy_landscape.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_dir / 'heteroplasmy_landscape.png'}")


def plot_hf_distribution(variants, out_dir):
    """Histogram of heteroplasmy frequencies."""
    low_level = [v["hf"] for v in variants if 0.001 < v["hf"] < 0.95]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.hist(low_level, bins=50, color="#4477AA", edgecolor="white", alpha=0.8)
    ax1.set_xlabel("Heteroplasmy Frequency")
    ax1.set_ylabel("Count")
    ax1.set_title("Distribution of Heteroplasmy Frequencies")
    ax1.axvline(0.01, color="red", linestyle="--", label="1% threshold")
    ax1.axvline(0.02, color="orange", linestyle="--", label="2% threshold")
    ax1.legend()

    ax2.hist(
        low_level,
        bins=np.logspace(np.log10(0.001), np.log10(0.95), 50),
        color="#44AA77",
        edgecolor="white",
        alpha=0.8,
    )
    ax2.set_xscale("log")
    ax2.set_xlabel("Heteroplasmy Frequency (log scale)")
    ax2.set_ylabel("Count")
    ax2.set_title("Distribution (Log Scale)")
    ax2.axvline(0.01, color="red", linestyle="--", label="1%")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(out_dir / "hf_distribution.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_dir / 'hf_distribution.png'}")


def generate_report(variants, sample_name, out_dir):
    """Generate markdown report on heteroplasmy findings."""
    low_level = [v for v in variants if v["hf"] < 0.95]
    homoplasmic = [v for v in variants if v["hf"] >= 0.95]

    age_hits = [v for v in low_level if v["pos"] in AGE_HOTSPOTS and v["hf"] > 0.003]
    pathogenic_hits = [
        v
        for v in variants
        if v["pos"] in PATHOGENIC
        and matches_pathogenic_change(v, PATHOGENIC[v["pos"]][0])
    ]
    above_1pct = [v for v in low_level if v["hf"] >= 0.01]
    above_2pct = [v for v in low_level if v["hf"] >= 0.02]
    hf_vals = [v["hf"] for v in low_level]

    lines = [
        "# Mitochondrial Heteroplasmy Report",
        "",
        f"**Sample:** {sample_name}",
        f"**Total variants called:** {len(variants)}",
        f"**Homoplasmic (HF \u2265 0.95):** {len(homoplasmic)}",
        f"**Heteroplasmic (HF < 0.95):** {len(low_level)}",
        f"**Mean read depth:** {np.mean([v['dp'] for v in variants]):.0f}x",
        "",
        "## Summary Statistics",
        "",
        f"- Variants above 1% HF: **{len(above_1pct)}**",
        f"- Variants above 2% HF: **{len(above_2pct)}**",
        f"- Median low-level HF: **{np.median(hf_vals):.4f}** ({np.median(hf_vals)*100:.2f}%)",
        f"- Max low-level HF: **{max(hf_vals):.4f}** at position m.{low_level[np.argmax(hf_vals)]['pos']}",
        "",
        "## Age-Associated Heteroplasmy Hotspots",
        "",
        "Several D-loop positions accumulate heteroplasmy with age. "
        "Here are the ones detected:",
        "",
        "| Position | Change | HF | Depth | Note |",
        "|----------|--------|-----|-------|------|",
    ]
    for v in sorted(age_hits, key=lambda x: -x["hf"]):
        lines.append(
            f"| m.{v['pos']} | {v['ref']}>{v['alt']} | {v['hf']:.3f} "
            f"| {v['dp']} | {AGE_HOTSPOTS[v['pos']]} |"
        )
    lines.extend(
        [
            "",
            "*Note: Low-level heteroplasmy at these positions is normal and increases "
            "with age. Values below ~5% are typical for younger individuals.*",
            "",
            "## MITOMAP Pathogenic Variant Screen",
            "",
        ]
    )

    if pathogenic_hits:
        lines.extend(
            [
                "\u26a0 **Variants found at known pathogenic positions:**",
                "",
                "| Position | Change | HF | Disease | Gene |",
                "|----------|--------|-----|---------|------|",
            ]
        )
        for v in pathogenic_hits:
            expected_change, disease, gene = PATHOGENIC[v["pos"]]
            lines.append(
                f"| m.{v['pos']} | {v['ref']}>{v['alt']} "
                f"(expected: {expected_change}) | {v['hf']:.4f} "
                f"| {disease} | {gene} |"
            )
        lines.extend(
            [
                "",
                "*Note: Very low heteroplasmy (<1%) at pathogenic positions is common "
                "in the general population and typically not clinically significant. "
                "These are sub-clonal variants.*",
            ]
        )
    else:
        lines.append(
            "No variants detected at confirmed MITOMAP pathogenic " "positions. \u2713"
        )
    lines.append("")

    lines.extend(
        [
            "## Top Heteroplasmic Variants (HF \u2265 1%)",
            "",
            "| Position | Change | HF | Depth | Region |",
            "|----------|--------|-----|-------|--------|",
        ]
    )
    for v in sorted(above_1pct, key=lambda x: -x["hf"]):
        name, _ = get_region_type(v["pos"])
        lines.append(
            f"| m.{v['pos']} | {v['ref']}>{v['alt']} "
            f"| {v['hf']:.3f} ({v['hf']*100:.1f}%) | {v['dp']} "
            f"| {name} |"
        )

    lines.extend(
        [
            "",
            "## Figures",
            "",
            "- `heteroplasmy_landscape.png` \u2014 Genome-wide heteroplasmy "
            "frequencies with gene map and depth",
            "- `hf_distribution.png` \u2014 Distribution of heteroplasmy frequencies",
            "",
            "---",
            "*Generated by mito_heteroplasmy.py*",
        ]
    )

    report = "\n".join(lines)
    (out_dir / "report.md").write_text(report)
    print(f"Saved: {out_dir / 'report.md'}")
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Mitochondrial heteroplasmy landscape analysis. "
        "Plots heteroplasmy frequencies, identifies age-related "
        "hotspots, and screens MITOMAP pathogenic mutations.",
    )
    parser.add_argument("vcf", help="Input mitochondrial VCF (may be gzipped)")
    parser.add_argument("output_dir", help="Directory for output report and plots")
    parser.add_argument(
        "--sample-name",
        default=None,
        help="Sample name for report header (default: read from VCF)",
    )
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_name = args.sample_name or get_sample_name_from_vcf(vcf_path)

    variants = extract_variants_from_vcf(vcf_path)

    plot_landscape(variants, out_dir)
    plot_hf_distribution(variants, out_dir)
    report = generate_report(variants, sample_name, out_dir)
    print()
    print(report)


if __name__ == "__main__":
    main()
