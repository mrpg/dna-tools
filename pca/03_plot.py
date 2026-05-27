#!/usr/bin/env python3
"""Plot ancestry PCA: user's genome vs. 1000 Genomes super-populations.

Usage: python3 dna-tools/pca/03_plot.py [outdir]

  [outdir]: directory containing pca.eigenvec and pca.eigenval
            (default: Reports/PCA)

Reads populations.tsv from dna-tools/pca/data/ (adjacent to this script).
Writes pca_pc1_pc2.png and pca_pc1_pc3.png to [outdir].
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPTDIR = Path(__file__).resolve().parent
REFDATA = SCRIPTDIR / "data"

SUPERPOP_META = {
    "AFR": ("African", "#E41A1C"),
    "AMR": ("Admixed American", "#377EB8"),
    "EAS": ("East Asian", "#4DAF4A"),
    "EUR": ("European", "#984EA3"),
    "SAS": ("South Asian", "#FF7F00"),
}


def load_populations():
    pops = {}
    with open(REFDATA / "populations.tsv") as f:
        next(f)  # header
        for line in f:
            sample, _pop, superpop, *_ = line.strip().split("\t")
            pops[sample] = superpop
    return pops


def load_eigenvec(results_dir):
    samples, pcs = [], []
    with open(results_dir / "pca.eigenvec") as f:
        next(f)  # header
        for line in f:
            parts = line.split()
            samples.append(parts[0])
            pcs.append([float(x) for x in parts[1:]])
    return samples, pcs


def load_eigenval(results_dir):
    return [
        float(line)
        for line in (results_dir / "pca.eigenval").read_text().splitlines()
        if line.strip()
    ]


def make_plot(samples, pcs, pops, pct_var, pc_x, pc_y, out_path):
    groups = {sp: ([], []) for sp in SUPERPOP_META}
    user_xy, user_name = None, None

    for sample, pc in zip(samples, pcs):
        if sample in pops and pops[sample] in groups:
            groups[pops[sample]][0].append(pc[pc_x])
            groups[pops[sample]][1].append(pc[pc_y])
        else:
            user_xy = (pc[pc_x], pc[pc_y])
            user_name = sample

    fig, ax = plt.subplots(figsize=(10, 8))
    for sp, (label, color) in SUPERPOP_META.items():
        if groups[sp][0]:
            ax.scatter(
                groups[sp][0],
                groups[sp][1],
                c=color,
                label=label,
                alpha=0.4,
                s=12,
                edgecolors="none",
            )
    if user_xy:
        ax.scatter(
            *user_xy,
            c="black",
            marker="*",
            s=350,
            zorder=10,
            label=f"You ({user_name})",
            edgecolors="white",
            linewidths=0.5,
        )

    ax.set_xlabel(f"PC{pc_x+1} ({pct_var[pc_x]:.1f}%)")
    ax.set_ylabel(f"PC{pc_y+1} ({pct_var[pc_y]:.1f}%)")
    ax.set_title("Ancestry PCA — 1000 Genomes + You")
    ax.legend(loc="best", framealpha=0.9, markerscale=2)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    plt.close(fig)

    return user_xy, user_name


def main():
    results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("Reports/PCA")

    if not (results_dir / "pca.eigenvec").exists():
        print(
            f"Error: {results_dir}/pca.eigenvec not found. Run 02_pca.sh first.",
            file=sys.stderr,
        )
        sys.exit(1)

    pops = load_populations()
    samples, pcs = load_eigenvec(results_dir)
    eigenvals = load_eigenval(results_dir)
    total = sum(eigenvals)
    pct_var = [100 * v / total for v in eigenvals]

    # PC1 vs PC2 (standard view)
    user_xy, user_name = make_plot(
        samples,
        pcs,
        pops,
        pct_var,
        0,
        1,
        results_dir / "pca_pc1_pc2.png",
    )

    # PC1 vs PC3 (often separates European sub-structure)
    make_plot(
        samples,
        pcs,
        pops,
        pct_var,
        0,
        2,
        results_dir / "pca_pc1_pc3.png",
    )

    if user_xy:
        print(f"\nYour sample: {user_name}")
        for i in range(min(5, len(pcs[0]))):
            for s, pc in zip(samples, pcs):
                if s == user_name:
                    print(f"  PC{i+1} = {pc[i]:.6f}")
                    break


if __name__ == "__main__":
    main()
