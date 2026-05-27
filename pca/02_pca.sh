#!/usr/bin/env bash
#
# 02_pca.sh — Ancestry PCA against 1000 Genomes
#
# Projects a single whole-genome sample onto a principal component
# analysis of 2504 unrelated individuals from the 1000 Genomes Phase 3
# panel (26 populations, 5 super-populations).
#
# The pipeline:
#   1. Converts 1000G per-chromosome VCFs to plink2 format, filtering
#      to common (MAF > 5%) biallelic SNPs.          [first run only]
#   2. Merges chromosomes into a single dataset.     [first run only]
#   3. Restricts to 2504 unrelated samples.          [first run only]
#   4. LD-prunes (window 50, step 5, r² < 0.2).     [first run only]
#   5. Extracts user genotypes from gVCF at the ~500k pruned sites,
#      correctly resolving GATK reference blocks to homozygous-ref.
#   6. Merges user + 1000G via bcftools and runs PCA (20 components).
#
# Steps 1–4 build shared reference data (stored in dna-tools/pca/data/)
# and only run once.  Steps 5–6 are per-sample.
#
# Designed for Sequencing.com 30x WGS output processed through gvcf.sh
# (GATK HaplotypeCaller gVCF, GRCh38, bare chromosome names 1–22/X/Y).
#
# Requirements:
#   - bcftools (tested with 1.21)
#   - bgzip, tabix (htslib)
#   - plink2  (tested with v2.0.0-a.6.9)
#   - Python 3 with pysam (tested with 0.23.3)
#   - 01_download.sh must have been run first
#
# Usage:
#   bash dna-tools/pca/02_pca.sh <gvcf> [outdir]
#
#   <gvcf>  : GATK HaplotypeCaller gVCF (GRCh38, .gz + .tbi indexed)
#   [outdir]: output directory (default: Reports/PCA)
#
# Output in [outdir]:
#   pca.eigenvec  — PC coordinates for all 2505 samples
#   pca.eigenval  — variance explained by each PC
#
# After the run, generate plots with:
#   python3 dna-tools/pca/03_plot.py [outdir]
#

set -euo pipefail

# ---- Argument parsing --------------------------------------------------------
if [ $# -lt 1 ]; then
    echo "Usage: $0 <gvcf> [outdir]" >&2
    exit 1
fi

USER_GVCF="$1"
OUTDIR="${2:-Reports/PCA}"

if [ ! -f "$USER_GVCF" ]; then
    echo "Error: gVCF not found: $USER_GVCF" >&2
    exit 1
fi

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
REFDATA="${SCRIPTDIR}/data"
KGDIR="${REFDATA}/1kgp"
PANEL="${REFDATA}/populations.tsv"
THREADS=4
MEM=4000  # MB for plink2

[[ -f "$PANEL" ]] || { echo "Missing $PANEL — run 01_download.sh first" >&2; exit 1; }

mkdir -p "$REFDATA" "$OUTDIR"

# ===========================================================
# Step 1: Convert each 1000G chromosome VCF → plink2
#         Keep biallelic SNPs with MAF > 5%, <1% missingness
# ===========================================================
echo "=== Step 1/6: Convert 1000G VCFs to plink2 ==="
MERGE_LIST="${REFDATA}/merge_list.txt"
> "$MERGE_LIST"

for chr in $(seq 1 22); do
    OUT="${REFDATA}/1kgp_chr${chr}"
    if [[ -f "${OUT}.pgen" ]]; then
        echo "${OUT}" >> "$MERGE_LIST"
        continue
    fi
    VCF="${KGDIR}/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    [[ -f "$VCF" ]] || { echo "Missing $VCF — run 01_download.sh first" >&2; exit 1; }
    echo "  chr${chr}..."
    plink2 --vcf "$VCF" \
        --snps-only just-acgt --max-alleles 2 --maf 0.05 --geno 0.01 \
        --set-all-var-ids '@:#:$r:$a' \
        --make-pgen --out "$OUT" \
        --threads "$THREADS" --memory "$MEM" 2>&1 | tail -1
    echo "${OUT}" >> "$MERGE_LIST"
done

# ===========================================================
# Step 2: Merge all 22 per-chromosome plink2 files
# ===========================================================
echo "=== Step 2/6: Merge chromosomes ==="
MERGED="${REFDATA}/1kgp_all"
if [[ ! -f "${MERGED}.pgen" ]]; then
    N_FILES=$(wc -l < "$MERGE_LIST")
    if (( N_FILES < 2 )); then
        echo "Need at least 2 chromosome files; found $N_FILES" >&2; exit 1
    fi
    plink2 --pmerge-list "$MERGE_LIST" pfile \
        --make-pgen --out "$MERGED" \
        --threads "$THREADS" --memory "$MEM" 2>&1 | tail -3
fi

# ===========================================================
# Step 3: Keep only the 2504 unrelated Phase 3 samples
# ===========================================================
echo "=== Step 3/6: Filter to 2504 unrelated samples ==="
FILTERED="${REFDATA}/1kgp_unrel"
if [[ ! -f "${FILTERED}.pgen" ]]; then
    awk 'NR > 1 { print $1 }' "$PANEL" > "${REFDATA}/keep_samples.txt"
    plink2 --pfile "$MERGED" \
        --keep "${REFDATA}/keep_samples.txt" \
        --make-pgen --out "$FILTERED" \
        --threads "$THREADS" 2>&1 | tail -1
fi

# ===========================================================
# Step 4: LD-prune (window 50, step 5, r² < 0.2)
# ===========================================================
echo "=== Step 4/6: LD pruning ==="
PRUNED="${REFDATA}/1kgp_pruned"
if [[ ! -f "${PRUNED}.pgen" ]]; then
    plink2 --pfile "$FILTERED" \
        --indep-pairwise 50 5 0.2 \
        --out "${REFDATA}/ld_prune" \
        --threads "$THREADS" 2>&1 | tail -1
    plink2 --pfile "$FILTERED" \
        --extract "${REFDATA}/ld_prune.prune.in" \
        --make-pgen --out "$PRUNED" \
        --threads "$THREADS" 2>&1 | tail -1
    echo "  Retained $(wc -l < "${REFDATA}/ld_prune.prune.in") SNPs after LD pruning"
fi

# ===========================================================
# Step 5: Extract user genotypes at pruned PCA sites
# ===========================================================
echo "=== Step 5/6: Extract user genotypes from gVCF ==="
USER_PLINK="${OUTDIR}/user"
if [[ ! -f "${USER_PLINK}.pgen" ]]; then
    echo "  Querying gVCF at $(grep -cv '^#' "${PRUNED}.pvar") sites..."
    python3 "${SCRIPTDIR}/extract_user_genotypes.py" \
        "$USER_GVCF" "${PRUNED}.pvar" "${OUTDIR}/user_pca.vcf"

    bgzip -c "${OUTDIR}/user_pca.vcf" > "${OUTDIR}/user_pca.vcf.gz"
    tabix -p vcf "${OUTDIR}/user_pca.vcf.gz"
    rm "${OUTDIR}/user_pca.vcf"

    plink2 --vcf "${OUTDIR}/user_pca.vcf.gz" \
        --snps-only just-acgt --max-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' \
        --make-pgen --out "$USER_PLINK" 2>&1 | tail -1
fi

# ===========================================================
# Step 6: Merge user + 1000G via bcftools, then PCA
#         (plink2 alpha's --pmerge doesn't support cross-sample merges yet)
# ===========================================================
echo "=== Step 6/6: Merge + PCA ==="
FINAL="${OUTDIR}/merged"
if [[ ! -f "${FINAL}.pgen" ]]; then
    echo "  Exporting 1000G pruned data to VCF..."
    plink2 --pfile "$PRUNED" \
        --export vcf bgz --out "${OUTDIR}/1kgp_pruned_export" \
        --threads "$THREADS" 2>&1 | tail -1
    tabix -p vcf "${OUTDIR}/1kgp_pruned_export.vcf.gz"

    echo "  Merging with user VCF..."
    bcftools merge \
        "${OUTDIR}/1kgp_pruned_export.vcf.gz" \
        "${OUTDIR}/user_pca.vcf.gz" \
        -Oz -o "${OUTDIR}/merged.vcf.gz"
    tabix -p vcf "${OUTDIR}/merged.vcf.gz"

    echo "  Converting to plink2..."
    plink2 --vcf "${OUTDIR}/merged.vcf.gz" \
        --make-pgen --out "$FINAL" \
        --threads "$THREADS" 2>&1 | tail -1

    # Remove intermediate VCFs (regenerated on re-run; inspect $OUTDIR before this line if debugging)
    rm -f "${OUTDIR}/1kgp_pruned_export.vcf.gz"* "${OUTDIR}/merged.vcf.gz"*
fi

plink2 --pfile "$FINAL" \
    --pca 20 \
    --out "${OUTDIR}/pca" \
    --threads "$THREADS" 2>&1 | tail -3

echo ""
echo "=== Done ==="
echo "  ${OUTDIR}/pca.eigenvec  — PC coordinates (2504 ref + 1 user)"
echo "  ${OUTDIR}/pca.eigenval  — variance explained per PC"
echo "Next: python3 dna-tools/pca/03_plot.py ${OUTDIR}"
