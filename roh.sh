#!/usr/bin/env bash
#
# roh.sh — Runs of Homozygosity (ROH) analysis
#
# Detects long homozygous stretches across autosomes to assess
# consanguinity or identical-by-descent (IBD) sharing. Two scans
# are performed:
#
#   1. Standard (>=500 kb) — segments at this length suggest recent
#      parental relatedness (consanguinity) when numerous/large.
#   2. Sensitive (>=300 kb) — captures shorter population-level IBD
#      segments common in all individuals.
#
# The inbreeding coefficient F(ROH) = total ROH length / autosomal
# genome length. Reference thresholds:
#
#   F(ROH) > 0.0625  →  parents are first cousins
#   F(ROH) > 0.0156  →  parents are second cousins
#   F(ROH) < 0.005   →  typical outbred European
#
# Designed for Sequencing.com 30x WGS output (GRCh38, VCFv4.2 with
# GT/DP/GQ/AD fields but no PL values in the SNP/indel genome VCF).
#
# Requirements:
#   - bcftools (tested with 1.21)
#   - plink2   (tested with v2.0.0-a.6.9; for VCF→BED conversion)
#   - plink1   (PLINK v1.07; has --homozyg, which plink2 lacks)
#
# Usage:
#   ./dna-tools/roh.sh <snp-indel-genome.vcf.gz> <sex> [outdir]
#
#   <sex>   : "male" or "female" (required for chrX handling)
#   [outdir]: output directory (default: Reports/ROH)
#
# Output:
#   roh_results.hom         — individual ROH segments (>=500 kb)
#   roh_results.hom.indiv   — per-sample summary (>=500 kb)
#   roh_sensitive.hom       — individual ROH segments (>=300 kb)
#   roh_sensitive.hom.indiv — per-sample summary (>=300 kb)
#   snps_only.vcf.gz        — intermediate: biallelic alt SNPs
#   snps_auto.{bed,bim,fam} — intermediate: autosomal PLINK fileset
#

set -euo pipefail

# ---- Argument parsing -----------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <snp-indel-genome.vcf.gz> <male|female> [outdir]" >&2
    exit 1
fi

VCF="$1"
SEX_STR="$2"
OUTDIR="${3:-Reports/ROH}"

if [ ! -f "$VCF" ]; then
    echo "Error: VCF not found: $VCF" >&2
    exit 1
fi

# PLINK sex codes: 1 = male, 2 = female
case "$SEX_STR" in
    male)   SEX_CODE=1 ;;
    female) SEX_CODE=2 ;;
    *)      echo "Error: sex must be 'male' or 'female'" >&2; exit 1 ;;
esac

# Auto-detect the sample name from the VCF header
SAMPLE="$(bcftools query -l "$VCF" | head -1)"
if [ -z "$SAMPLE" ]; then
    echo "Error: no samples found in VCF" >&2
    exit 1
fi
echo "Detected sample: $SAMPLE (sex: $SEX_STR)"

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Step 1: Extract biallelic SNPs that carry at least one alt allele.
#
# Sequencing.com genome VCFs contain ~12M sites, most homozygous-reference.
# bcftools roh segfaults on these (missing PL fields + --GTs-only bug in
# v1.21), so we filter down to ~4M informative biallelic SNPs for PLINK.
# ---------------------------------------------------------------------------
echo "Step 1: Extracting biallelic alt SNPs from genome VCF..."
bcftools view -v snps -m2 -M2 -i 'GT="alt"' "$VCF" \
    -Oz -o "$OUTDIR/snps_only.vcf.gz"
bcftools index "$OUTDIR/snps_only.vcf.gz"

# ---------------------------------------------------------------------------
# Step 2: Convert VCF to PLINK binary format.
#
# --update-sex: required because chrX is present and plink2 refuses to
#   process it without sex information.
# --split-par hg38: reclassifies pseudoautosomal region SNPs on chrX so
#   they are handled correctly (diploid, not haploid).
# --allow-extra-chr: the VCF may contain non-standard contigs (MT,
#   scaffolds) that would otherwise cause an error; dropped in Step 3.
# ---------------------------------------------------------------------------
echo "Step 2: Converting to PLINK binary format..."
echo "0 $SAMPLE $SEX_CODE" > "$OUTDIR/sex.txt"
plink2 --vcf "$OUTDIR/snps_only.vcf.gz" \
    --make-bed \
    --out "$OUTDIR/snps_plink" \
    --allow-extra-chr \
    --update-sex "$OUTDIR/sex.txt" \
    --split-par hg38

# ---------------------------------------------------------------------------
# Step 3: Restrict to autosomes (chr 1-22).
#
# PLINK v1.07 expects numeric chromosome codes and cannot parse "X", "MT",
# "PAR1", etc. ROH analysis on sex chromosomes is also not meaningful for
# assessing consanguinity, so we simply exclude them.
# ---------------------------------------------------------------------------
echo "Step 3: Filtering to autosomes (chr 1-22)..."
plink2 --bfile "$OUTDIR/snps_plink" \
    --chr 1-22 \
    --make-bed \
    --out "$OUTDIR/snps_auto"

# ---------------------------------------------------------------------------
# Step 4a: Standard ROH scan (>=500 kb).
#
# Parameters follow Howrigan et al. (2011) / Pemberton et al. (2012)
# conventions for WGS-density data:
#   --homozyg-snp 50    : require >=50 SNPs per segment
#   --homozyg-kb 500    : minimum segment length 500 kb
#   --homozyg-density 50: max 50 kb per SNP (ensures adequate marker density)
#   --homozyg-gap 1000  : max 1 Mb gap between consecutive SNPs
#   --homozyg-window-snp 50  : sliding window size
#   --homozyg-window-het 1   : allow 1 het per window (genotyping errors)
#   --homozyg-window-missing 5: allow 5 missing calls per window
# ---------------------------------------------------------------------------
echo "Step 4a: ROH detection (>=500 kb, consanguinity-level)..."
plink1 --bfile "$OUTDIR/snps_auto" \
    --homozyg \
    --homozyg-snp 50 \
    --homozyg-kb 500 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --homozyg-window-snp 50 \
    --homozyg-window-het 1 \
    --homozyg-window-missing 5 \
    --out "$OUTDIR/roh_results" \
    --noweb

# ---------------------------------------------------------------------------
# Step 4b: Sensitive ROH scan (>=300 kb).
#
# Lower thresholds capture shorter IBD segments from more distant shared
# ancestry. These are normal in all individuals and useful for comparing
# population-level homozygosity patterns.
# ---------------------------------------------------------------------------
echo "Step 4b: ROH detection (>=300 kb, population-level IBD)..."
plink1 --bfile "$OUTDIR/snps_auto" \
    --homozyg \
    --homozyg-snp 30 \
    --homozyg-kb 300 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --homozyg-window-snp 30 \
    --homozyg-window-het 1 \
    --homozyg-window-missing 5 \
    --out "$OUTDIR/roh_sensitive" \
    --noweb

# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------
echo ""
echo "========================================="
echo "  ROH Analysis Complete"
echo "========================================="
echo ""
echo "Standard scan (>=500 kb):"
cat "$OUTDIR/roh_results.hom.indiv"
echo ""
echo "Sensitive scan (>=300 kb):"
cat "$OUTDIR/roh_sensitive.hom.indiv"
echo ""
echo "Detailed segment lists: $OUTDIR/roh_results.hom, $OUTDIR/roh_sensitive.hom"
