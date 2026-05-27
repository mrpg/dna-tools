#!/usr/bin/env python3
"""Extract user genotypes from a gVCF at specific PCA sites.

Reads a plink2 .pvar file for the target positions, then queries
the gVCF at each one.  Reference blocks are correctly resolved to
0/0 (homozygous reference).

Usage: python3 extract_user_genotypes.py <gvcf> <pvar> <output.vcf>
"""

import sys
from time import time

import pysam


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <gvcf> <pvar> <output.vcf>", file=sys.stderr)
        sys.exit(1)

    gvcf_path, pvar_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    # Read target sites from plink2 .pvar (tab-separated, #-header lines)
    sites = []
    with open(pvar_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt = line.strip().split("\t")[:5]
            sites.append((chrom, int(pos), vid, ref, alt))

    vcf_in = pysam.VariantFile(gvcf_path)
    sample_name = list(vcf_in.header.samples)[0]

    ref_ref = variant = missing = 0
    total = len(sites)
    t0 = time()

    with open(out_path, "w") as out:
        # Minimal VCF header
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for i in range(1, 23):
            out.write(f"##contig=<ID={i}>\n")
        out.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )

        for idx, (chrom, pos, vid, ref_al, alt_al) in enumerate(sites):
            if idx % 20000 == 0 and idx > 0:
                elapsed = time() - t0
                eta = elapsed / idx * (total - idx)
                print(
                    f"  {idx}/{total} ({100*idx/total:.0f}%) — ~{eta:.0f}s left",
                    flush=True,
                )

            gt_str = "./."

            try:
                records = list(vcf_in.fetch(str(chrom), pos - 1, pos))
            except (ValueError, KeyError):
                records = []

            hit = False
            for rec in records:
                rec_start = rec.pos  # 1-based VCF POS (pysam 0.23+)
                rec_end = rec.stop  # 1-based inclusive end

                if not (rec_start <= pos <= rec_end):
                    continue

                gt = rec.samples[0]["GT"]
                if gt is None or all(a is None for a in gt):
                    hit = True
                    missing += 1
                    break

                if rec_start == pos:
                    # Record at exact target position — could be a real variant
                    real_alts = [
                        str(a) for a in rec.alleles[1:] if a and str(a) != "<NON_REF>"
                    ]
                    if alt_al in real_alts:
                        alt_idx = list(rec.alleles).index(alt_al)
                        gt_str = "/".join("1" if a == alt_idx else "0" for a in gt)
                        variant += 1
                    else:
                        # Different alt or ref-only — user is ref/ref for this allele
                        gt_str = "0/0"
                        ref_ref += 1
                else:
                    # Reference block encompassing the target position
                    gt_str = "0/0"
                    ref_ref += 1

                hit = True
                break

            if not hit:
                missing += 1

            out.write(
                f"{chrom}\t{pos}\t{vid}\t{ref_al}\t{alt_al}\t.\t.\t.\tGT\t{gt_str}\n"
            )

    vcf_in.close()

    print(
        f"Done — {total} sites: {variant} variant, {ref_ref} ref/ref, "
        f"{missing} missing ({100*missing/total:.1f}%)"
    )


if __name__ == "__main__":
    main()
