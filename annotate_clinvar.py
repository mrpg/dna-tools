#!/usr/bin/env python3
"""
Annotate personal gVCF variants against ClinVar to find clinically relevant variants.

Usage:
    python annotate_clinvar.py <sample.g.vcf.gz> <clinvar.vcf.gz> [output.tsv] [reference.fa] > report.md

The sample file must be a gVCF (genomic VCF) containing multi-allelic ALT
records — standard single-ALT VCFs are not supported because this script
relies on multi-allelic sites to correctly resolve genotype allele indices.

Both files must use the same chromosome naming convention (e.g., both
Ensembl-style "1, 2, ..." or both UCSC-style "chr1, chr2, ...") and the
same reference build (GRCh38).

The clinvar.vcf.gz file must be accompanied by a clinvar.vcf.gz.tbi file.
    tabix -p vcf clinvar.vcf.gz

Matching is by exact (chrom, pos, ref, alt), and the sample's genotype is
checked to confirm the matched ALT allele is actually carried. Each sample
allele is normalized before lookup: multi-allelic gVCF records pad every
allele to the longest REF (e.g. TAAA>TAA for a 1-bp deletion), so the shared
suffix/prefix is trimmed and, if reference.fa is given, indels are
left-aligned. Output pos/ref/alt use the normalized allele; decoded_genotype
keeps the record's original allele spelling. Always pass reference.fa (the
same FASTA used for ClinVar) for complete indel matching. ClinVar must be
normalized beforehand:
    bcftools norm \
            -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT \
            -f Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            -o clinvar_norm.vcf.gz -O z clinvar.vcf.gz

After that, reindex with:
    tabix -p vcf clinvar_norm.vcf.gz
"""

import sys

import pandas as pd
import pysam
from tqdm import tqdm


def count_vcf_records(path):
    """Count records in a VCF file for progress bar sizing."""
    vcf = pysam.VariantFile(path)
    n = 0
    for _ in vcf:
        n += 1
    vcf.close()
    return n


def parse_clinvar(clinvar_path):
    """Build lookup from ClinVar: (chrom, pos, ref, alt) -> info dict."""
    total = count_vcf_records(clinvar_path)
    clinvar = {}
    vcf = pysam.VariantFile(clinvar_path)
    for rec in tqdm(vcf, total=total, desc="Loading ClinVar", unit=" records"):
        for alt in rec.alts or []:
            key = (rec.chrom, rec.pos, rec.ref, alt)
            info = rec.info
            clnsig = info.get("CLNSIG", "")
            if isinstance(clnsig, tuple):
                clnsig = ",".join(clnsig)
            clinvar[key] = {
                "clnsig": clnsig,
                "clndn": info.get("CLNDN", ""),
                "geneinfo": info.get("GENEINFO", ""),
                "clnrevstat": info.get("CLNREVSTAT", ""),
                "clnhgvs": info.get("CLNHGVS", ""),
                "rs": info.get("RS", ""),
                "mc": info.get("MC", ""),
            }
            for k, v in clinvar[key].items():
                if isinstance(v, tuple):
                    clinvar[key][k] = ",".join(str(x) for x in v)
    print(f"Loaded {len(clinvar):,} ClinVar variant entries.", file=sys.stderr)
    return clinvar


def normalize_allele(chrom, pos, ref, alt, fasta=None):
    """Trim shared bases and left-align (vt-style) one REF/ALT pair.

    Without a FASTA, or on a contig the FASTA lacks, indels are only trimmed.
    """
    orig = (pos, ref, alt)
    ref, alt = ref.upper(), alt.upper()
    changed = True
    while changed:
        changed = False
        if ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            if fasta is None and (len(ref) == 1 or len(alt) == 1):
                break
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if (not ref or not alt) and pos > 1:
            try:
                base = fasta.fetch(chrom, pos - 2, pos - 1).upper()
            except KeyError:
                return normalize_allele(chrom, *orig)
            ref, alt, pos = base + ref, base + alt, pos - 1
            changed = True
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def annotate_sample(sample_path, clinvar, fasta=None):
    """Scan sample VCF and match against ClinVar."""
    total = count_vcf_records(sample_path)
    vcf = pysam.VariantFile(sample_path)
    results = []
    matched = 0
    renormalized = 0

    for rec in tqdm(vcf, total=total, desc="Scanning sample", unit=" variants"):
        if rec.alts is None or rec.alts[0] == ".":
            continue
        if rec.filter.keys() and "PASS" not in rec.filter.keys():
            continue

        sample = rec.samples[0]
        gt = sample.get("GT", (None, None))
        dp = sample.get("DP", None)
        gq = sample.get("GQ", None)

        for alt_idx, alt in enumerate(rec.alts, start=1):
            if alt in (".", "*", "<NON_REF>", "<*>"):
                continue
            # Only carried alleles can match; skip the normalization cost otherwise
            if gt and alt_idx not in gt:
                continue
            pos, ref, norm_alt = normalize_allele(
                rec.chrom, rec.pos, rec.ref, alt, fasta
            )
            key = (rec.chrom, pos, ref, norm_alt)
            if key in clinvar:
                matched += 1
                if (pos, ref, norm_alt) != (rec.pos, rec.ref, alt):
                    renormalized += 1
                cv = clinvar[key]
                gt_str = (
                    "/".join(str(a) if a is not None else "." for a in gt)
                    if gt
                    else "."
                )

                gene = cv["geneinfo"].split(":")[0] if cv["geneinfo"] else ""

                mc = cv["mc"]
                if mc:
                    mc = ",".join(
                        p.split("|")[-1] if "|" in p else p for p in mc.split(",")
                    )

                revstat = cv["clnrevstat"]
                stars = 0
                if "practice_guideline" in revstat:
                    stars = 4
                elif "reviewed_by_expert_panel" in revstat:
                    stars = 3
                elif "multiple_submitters" in revstat:
                    stars = 2
                elif "single_submitter" in revstat:
                    stars = 1

                # Build allele list for readable genotype display
                allele_map = [rec.ref] + list(rec.alts or [])
                decoded_gt = (
                    "/".join(
                        allele_map[a] if a is not None and a < len(allele_map) else "?"
                        for a in gt
                    )
                    if gt
                    else "."
                )

                results.append(
                    {
                        "chrom": rec.chrom,
                        "pos": pos,
                        "rsid": (
                            rec.id if rec.id else (f"rs{cv['rs']}" if cv["rs"] else "")
                        ),
                        "ref": ref,
                        "alt": norm_alt,
                        "genotype": gt_str,
                        "decoded_genotype": decoded_gt,
                        "depth": dp,
                        "qual": gq,
                        "gene": gene,
                        "significance": cv["clnsig"],
                        "condition": cv["clndn"].replace("_", " "),
                        "consequence": mc,
                        "hgvs": cv["clnhgvs"],
                        "review_stars": stars,
                        "review_status": revstat,
                    }
                )

    print(
        f"Done: {matched:,} ClinVar matches "
        f"({renormalized:,} only after allele normalization).",
        file=sys.stderr,
    )
    return results


def format_genotype(row):
    """Return a human-readable genotype description."""
    gt = row["genotype"]
    decoded_gt = row.get("decoded_genotype", "")
    alleles = gt.replace("|", "/").split("/")

    if decoded_gt:
        genotype_str = decoded_gt
    else:
        ref, alt = row["ref"], row["alt"]
        decoded = []
        for a in alleles:
            if a == "0":
                decoded.append(ref)
            elif a == "1":
                decoded.append(alt)
            elif a == ".":
                decoded.append("?")
            else:
                decoded.append(a)
        genotype_str = decoded_gt = "/".join(decoded)

    # Determine zygosity from the decoded alleles
    decoded_alleles = genotype_str.split("/")
    unique = set(decoded_alleles)
    if len(unique) == 1 and "?" not in unique:
        zygosity = "homozygous"
    elif len(unique) > 1:
        zygosity = "heterozygous"
    else:
        zygosity = None
    if zygosity:
        return f"{genotype_str} ({zygosity})"
    return genotype_str


def format_stars(n):
    """Return a star rating string like '\\*\\*--' (escaped for markdown)."""
    return r"\*" * n + "-" * (4 - n)


def print_variant_block(row, show_depth=False):
    """Print a single variant as a markdown list item."""
    depth_str = (
        f", depth={row['depth']}" if show_depth and pd.notna(row["depth"]) else ""
    )
    rsid_str = f" — {row['rsid']}" if row["rsid"] else ""

    print(
        f"- **{row['gene'] or '?'}** [{format_stars(row['review_stars'])}]"
        f"{rsid_str}"
    )
    print(f"  Genotype: {format_genotype(row)}{depth_str}")
    print(f"  {row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}")
    print(f"  *{row['condition']}*")
    if row["consequence"]:
        print(f"  {row['consequence']}")
    print()


def print_report(df):
    """Print a Pandoc-compatible markdown report to stdout."""
    print("# ClinVar Annotation Report")
    print()

    # Summary table
    print("## Summary by Clinical Significance")
    print()
    counts = df["significance"].value_counts()
    print("| Significance | Count |")
    print("|:-------------|------:|")
    for sig, count in counts.items():
        print(f"| {sig} | {count:,} |")
    print()

    # One section per significance category, ordered by max stars then count
    sig_order = (
        df.groupby("significance")["review_stars"]
        .max()
        .sort_values(ascending=False)
        .index
    )
    for sig in sig_order:
        subset = df[df["significance"] == sig].sort_values(
            ["review_stars", "gene"], ascending=[False, True]
        )
        count = len(subset)
        label = sig.replace("_", " ").replace("/", " / ").title() if sig else "Unknown"
        print(f"## {label} ({count:,})")
        print()
        for _, row in subset.iterrows():
            print_variant_block(row, show_depth=True)


def main():
    if len(sys.argv) < 3:
        print(__doc__.strip(), file=sys.stderr)
        sys.exit(1)

    sample_path = sys.argv[1]
    clinvar_path = sys.argv[2]
    output_path = sys.argv[3] if len(sys.argv) > 3 else "clinvar_annotated_variants.tsv"
    fasta = pysam.FastaFile(sys.argv[4]) if len(sys.argv) > 4 else None
    if fasta is not None:
        with pysam.VariantFile(sample_path) as vcf:
            if not set(fasta.references) & set(vcf.header.contigs):
                sys.exit(
                    "Error: reference.fa shares no contig names with the sample VCF "
                    "(e.g. 'chr1' vs '1')."
                )
    else:
        print(
            "Warning: no reference.fa given; indels are trimmed but not left-aligned.",
            file=sys.stderr,
        )

    clinvar = parse_clinvar(clinvar_path)
    results = annotate_sample(sample_path, clinvar, fasta)

    df = pd.DataFrame(results)
    if df.empty:
        print("No ClinVar matches found!", file=sys.stderr)
        return

    df.to_csv(output_path, sep="\t", index=False)
    print(
        f"Full results saved to {output_path} ({len(df):,} variants).", file=sys.stderr
    )

    print_report(df)


if __name__ == "__main__":
    main()
