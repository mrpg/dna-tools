#!/usr/bin/env bash
#
# 01_download.sh — Download 1000 Genomes Phase 3 reference data for PCA
#
# Downloads 22 autosomal VCFs (GRCh38, 2548 samples, biallelic SNPs +
# indels) from the International Genome Sample Resource (IGSR), plus the
# population panel mapping each sample to one of 26 populations and 5
# super-populations (AFR, AMR, EAS, EUR, SAS).
#
# Total download: ~12.7 GB.  Uses wget -c, so re-running safely resumes
# interrupted downloads.  These files are shared reference data — they
# only need to be downloaded once and are reused across samples.
#
# Source:
#   VCFs:  ftp.1000genomes.ebi.ac.uk/.../20190312_biallelic_SNV_and_INDEL/
#   Panel: ftp.1000genomes.ebi.ac.uk/.../integrated_call_samples_v3.20130502.ALL.panel
#
# Requirements:
#   - wget
#
# Usage:
#   bash dna-tools/pca/01_download.sh
#
# Output:
#   dna-tools/pca/data/1kgp/ALL.chr*.vcf.gz{,.tbi}  — per-chromosome VCFs
#   dna-tools/pca/data/populations.tsv                — sample→population map
#

set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
KGDIR="${DIR}/data/1kgp"
mkdir -p "$KGDIR"

BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"

echo "=== Downloading 1000 Genomes Phase 3 (GRCh38) — ~12.7 GB total ==="
for chr in $(seq 1 22); do
    F="ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    if [[ -f "${KGDIR}/${F}" && -f "${KGDIR}/${F}.tbi" ]]; then
        echo "chr${chr}: already downloaded"
        continue
    fi
    echo "--- chr${chr} ---"
    wget -c -q --show-progress -P "$KGDIR" "${BASE}/${F}"
    wget -c -q --show-progress -P "$KGDIR" "${BASE}/${F}.tbi"
done

echo ""
echo "=== Downloading population metadata ==="
PANEL="${DIR}/data/populations.tsv"
if [[ ! -f "$PANEL" ]]; then
    wget -q -O "$PANEL" \
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    echo "Saved: $PANEL"
else
    echo "Already exists: $PANEL"
fi

echo ""
echo "=== Download complete ==="
echo "Next: bash dna-tools/pca/02_pca.sh <your-gvcf>"
echo "After 02_pca.sh finishes, you can delete dna-tools/pca/data/1kgp/ to reclaim ~12.7 GB."
