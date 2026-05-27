#!/usr/bin/env python3
"""
Mitochondrial haplogroup screening from heteroplasmy-aware VCF.

Uses a small local subset of PhyloTree Build 17 defining mutations to identify
the best-matching candidate from the haplogroups listed in this file.

Usage:
    python mito_haplogroup.py <mito.vcf.gz> <output_dir> [--sample-name NAME]

Input VCF must have FORMAT fields: GT, HF (heteroplasmy frequency).
"""

import argparse
import subprocess
import sys
from pathlib import Path

# PhyloTree Build 17 haplogroup-defining mutations (rCRS-relative)
# Each entry: haplogroup -> set of (pos, alt) defining that branch
PHYLOTREE = {
    # Macro-haplogroup N (out of Africa, non-M branch)
    "L3": {(769, "G"), (1018, "G"), (16311, "C")},
    "N": {(8701, "G"), (9540, "C"), (10398, "G"), (10873, "C"), (15301, "A")},
    "R": {(12705, "T"), (16223, "T")},
    # European haplogroups
    "HV": {(11719, "A"), (14766, "T")},
    "H": {(7028, "C")},
    "H1": {(3010, "A")},
    "H1e": {(8473, "C")},
    "H1e1": {(16311, "C")},
    "H1e1a": {(477, "C"), (16519, "C")},
    "H2": {(1438, "G"), (4769, "G")},
    "H5": {(456, "T"), (16304, "C")},
    "H3": {(6776, "C")},
    "H7": {(4793, "G")},
    "V": {(4580, "A"), (15904, "T")},
    "HV0": {(16298, "C")},
    # Other major European
    "J": {(295, "T"), (489, "C"), (10398, "G"), (13708, "A"), (16069, "T")},
    "J1": {(3010, "A"), (462, "T")},
    "T": {
        (709, "A"),
        (1888, "A"),
        (4917, "G"),
        (10463, "C"),
        (13368, "A"),
        (14905, "A"),
        (15607, "G"),
        (15928, "A"),
    },
    "T2": {(11812, "G"), (14233, "G")},
    "U": {(11467, "G"), (12308, "G"), (12372, "A")},
    "U5": {(3197, "C"), (9477, "A"), (13617, "C")},
    "K": {(1189, "C"), (1811, "G"), (9055, "A"), (10550, "G")},
    "I": {(199, "C"), (204, "C"), (250, "C"), (10034, "C"), (10238, "C"), (16129, "A")},
    "W": {
        (119, "C"),
        (195, "C"),
        (204, "C"),
        (207, "A"),
        (1243, "C"),
        (3505, "G"),
        (5460, "A"),
        (8251, "A"),
        (8994, "A"),
        (11947, "G"),
        (15884, "C"),
    },
    "X": {(1719, "A"), (6221, "C"), (6371, "T"), (13966, "G"), (14470, "C")},
}

# Candidate branches in this limited local subset.
TREE = [
    ("H1e1a", ["H", "H1", "H1e", "H1e1", "H1e1a"]),
    ("H1e1", ["H", "H1", "H1e", "H1e1"]),
    ("H1e", ["H", "H1", "H1e"]),
    ("H1", ["H", "H1"]),
    ("H2", ["H", "H2"]),
    ("H3", ["H", "H3"]),
    ("H5", ["H", "H5"]),
    ("H7", ["H", "H7"]),
    ("H", ["H"]),
    ("V", ["HV", "HV0", "V"]),
    ("J1", ["J", "J1"]),
    ("J", ["J"]),
    ("T2", ["T", "T2"]),
    ("T", ["T"]),
    ("U5", ["U", "U5"]),
    ("K", ["U", "K"]),
    ("U", ["U"]),
    ("I", ["I"]),
    ("W", ["W"]),
    ("X", ["X"]),
]

# Mitochondrial genome regions (rCRS annotation)
MITO_REGIONS = [
    (1, 576, "D-loop / Control Region"),
    (577, 647, "tRNA-Phe"),
    (648, 1601, "12S rRNA"),
    (1602, 1670, "tRNA-Val"),
    (1671, 3229, "16S rRNA"),
    (3230, 3304, "tRNA-Leu(UUR)"),
    (3307, 4262, "ND1"),
    (4263, 4331, "tRNA-Ile"),
    (4329, 4400, "tRNA-Gln"),
    (4402, 4469, "tRNA-Met"),
    (4470, 5511, "ND2"),
    (5512, 5579, "tRNA-Trp"),
    (5587, 5655, "tRNA-Ala"),
    (5657, 5729, "tRNA-Asn"),
    (5761, 5826, "tRNA-Cys"),
    (5826, 5891, "tRNA-Tyr"),
    (5904, 7445, "COX1"),
    (7446, 7514, "tRNA-Ser(UCN)"),
    (7518, 7585, "tRNA-Asp"),
    (7586, 8269, "COX2"),
    (8295, 8364, "tRNA-Lys"),
    (8366, 8572, "ATP8"),
    (8527, 9207, "ATP6"),
    (9207, 9990, "COX3"),
    (9991, 10058, "tRNA-Gly"),
    (10059, 10404, "ND3"),
    (10405, 10469, "tRNA-Arg"),
    (10470, 10766, "ND4L"),
    (10760, 12137, "ND4"),
    (12138, 12206, "tRNA-His"),
    (12207, 12265, "tRNA-Ser(AGY)"),
    (12266, 12336, "tRNA-Leu(CUN)"),
    (12337, 14148, "ND5"),
    (14149, 14673, "ND6"),
    (14674, 14742, "tRNA-Glu"),
    (14747, 15887, "CytB"),
    (15888, 15953, "tRNA-Thr"),
    (15956, 16023, "tRNA-Pro"),
    (16024, 16569, "D-loop / Control Region"),
]


def get_region(pos):
    for start, end, name in MITO_REGIONS:
        if start <= pos <= end:
            return name
    return "Intergenic"


def extract_variants_from_vcf(vcf_path):
    """Extract variants from a mito VCF using bcftools. Returns pos -> [(alt, hf)]."""
    result = subprocess.run(
        ["bcftools", "query", "-f", r"%POS\t%REF\t%ALT\t[%HF]\n", str(vcf_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    variants = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        pos = int(parts[0])
        alts = parts[2].split(",")
        hfs = parts[3].split(",")
        for alt, hf in zip(alts, hfs):
            variants.setdefault(pos, []).append((alt, float(hf)))
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


def check_mutation(variants, pos, alt, threshold=0.80):
    """Check if a specific mutation is present above threshold."""
    if pos in variants:
        for v_alt, hf in variants[pos]:
            if v_alt == alt and hf >= threshold:
                return True, hf
    return False, 0.0


def determine_haplogroup(variants):
    """Walk the PhyloTree and find the best-matching haplogroup."""
    results = []
    for hg_name, path in TREE:
        all_mutations = set()
        for step in path:
            if step in PHYLOTREE:
                all_mutations.update(PHYLOTREE[step])

        matched, missed = [], []
        for pos, alt in sorted(all_mutations):
            found, hf = check_mutation(variants, pos, alt)
            if found:
                matched.append((pos, alt, hf))
            else:
                missed.append((pos, alt))

        if all_mutations:
            score = len(matched) / len(all_mutations)
            results.append((hg_name, score, matched, missed, path))

    results.sort(key=lambda x: (x[1], len(x[4])), reverse=True)
    return results


def build_report(variants, sample_name):
    """Build the haplogroup report as a string."""
    homoplasmic = []
    for pos, alts in sorted(variants.items()):
        for alt, hf in alts:
            if hf >= 0.95:
                homoplasmic.append((pos, alt, hf))

    results = determine_haplogroup(variants)
    best = results[0]

    lines = [
        "# Mitochondrial Haplogroup Report",
        "",
        f"**Sample:** {sample_name}",
        "**Reference:** rCRS (revised Cambridge Reference Sequence)",
        "**Method:** Limited local subset of PhyloTree Build 17 defining mutations",
        "**Scope:** Candidate screen only; confirm with a full mitochondrial haplogroup tool before using as a final assignment.",
        "",
        "## Result",
        "",
        f"**Best-Matching Candidate: {best[0]}**",
        f"- Match score: {best[1]:.0%} ({len(best[2])}/{len(best[2])+len(best[3])} defining mutations found)",
        f"- Candidate branch path in local subset: {' \u2192 '.join(best[4])}",
        "",
        "## Haplogroup Context",
        "",
    ]

    if best[0].startswith("H1e"):
        lines.extend(
            [
                "Haplogroup **H1e** is a subclade of H1, itself the most common subclade of "
                "haplogroup H \u2014 the dominant mitochondrial lineage in Europe (~40-50% of Europeans). "
                "H1 is particularly prevalent in Western Europe, with highest frequencies in Iberia "
                "and southwestern France.",
                "",
                "H1e is defined by the mutation T8473C and is found across Western and Central Europe. "
                "The H1e1 subclade (with T16311C) and H1e1a (adding T477C and T16519C) represent "
                "further refinements within this lineage.",
                "",
                "The H1 branch is estimated to have expanded ~10,000-13,000 years ago during the "
                "post-glacial recolonization of Europe from Franco-Cantabrian refugia.",
            ]
        )
    else:
        lines.append(
            f"Candidate {best[0]} \u2014 see PhyloTree (https://www.phylotree.org/) or a full haplogroup caller for detailed context."
        )
    lines.append("")

    lines.extend(
        [
            "## Homoplasmic Variants (HF \u2265 0.95)",
            "",
            "Fixed mitochondrial differences from the rCRS reference:",
            "",
            "| Position | Change | Heteroplasmy | Region |",
            "|----------|--------|-------------|--------|",
        ]
    )
    for pos, alt, hf in homoplasmic:
        lines.append(f"| m.{pos} | {alt} | {hf:.3f} | {get_region(pos)} |")

    lines.extend(
        [
            "",
            "## Defining Mutation Match Detail",
            "",
            f"### {best[0]} \u2014 matched mutations",
            "",
        ]
    )
    for pos, alt, hf in sorted(best[2]):
        lines.append(f"- m.{pos}{alt} (HF={hf:.3f}) \u2713")
    if best[3]:
        lines.extend(["", f"### {best[0]} \u2014 missing mutations", ""])
        for pos, alt in sorted(best[3]):
            lines.append(f"- m.{pos}{alt} \u2717")

    lines.extend(
        [
            "",
            "## Top Haplogroup Candidates",
            "",
            "| Haplogroup | Score | Matched | Missed |",
            "|------------|-------|---------|--------|",
        ]
    )
    for hg, score, matched, missed, path in results[:10]:
        lines.append(f"| {hg} | {score:.0%} | {len(matched)} | {len(missed)} |")

    lines.extend(
        [
            "",
            "---",
            "*Generated by mito_haplogroup.py \u2014 limited local matcher*",
        ]
    )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Determine mitochondrial haplogroup from a heteroplasmy-aware VCF.",
    )
    parser.add_argument("vcf", help="Input mitochondrial VCF (may be gzipped)")
    parser.add_argument("output_dir", help="Directory for output report")
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
    report = build_report(variants, sample_name)

    out_path = out_dir / "report.md"
    out_path.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
