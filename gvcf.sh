#!/usr/bin/env bash

set -euo pipefail

# It is recommended to run this script interactively (that is, step by step).
# If you are certain to not do that, remove this:

echo "Please run this script MANUALLY, step by step." > /dev/stderr && exit 41

# Generate gVCF from BAM on Debian/Ubuntu
# Requires: BAM file, reference genome (GRCh38)

# Walkthrough placeholder: replace this with your *.bam path before running.
BAMFILE="your-sample.bam"

[ -f "${BAMFILE}" ] || {
    echo "Error. Please replace BAMFILE with your BAM path in this walkthrough script." > /dev/stderr
    exit 42
}

# 1. Install dependencies
sudo apt update
sudo apt install -y openjdk-25-jre samtools

# 2. Download GATK (check https://github.com/broadinstitute/gatk/releases for latest)
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip
export PATH="$PWD/gatk-4.6.2.0:$PATH"

# 3. Download GRCh38 reference (if you don't have it)
# This is ~841 MiB compressed, ~2.9 GiB uncompressed

wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 4. Index the reference (required by GATK)
# This takes few seconds
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.primary_assembly.fa

# 5. Check your BAM file
# Should output OK
samtools quickcheck "${BAMFILE}" && echo OK

# 6. Index your BAM file
# This takes about 2 minutes
samtools index "${BAMFILE}"

# 7. Generate gVCF using HaplotypeCaller
# This takes about 10 hours
# You can adjust the number for -native-pair-hmm-threads to change
# the degree of parallelism
gatk HaplotypeCaller \
    -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -I "${BAMFILE}" \
    -O output.g.vcf.gz \
    -ERC GVCF \
    -native-pair-hmm-threads 8

# The output.g.vcf.gz file is your gVCF!
