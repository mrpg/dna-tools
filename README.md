# dna-tools

Simple scripts to work with 30x whole genome sequencing data (GRCh38).

**Do not run these scripts blindly.** Read each one, understand what it does, and adapt paths and parameters to your own data before executing. These are working notes shared as-is — not a turnkey pipeline. You are responsible for verifying every result independently.

Tested with data from [Sequencing.com](https://sequencing.com/). Other providers may or may not work. There is no warranty. Not affiliated with Sequencing.com.

[Read more here.](https://max.pm/posts/sequencing/)

## What's included

| Script | Description |
|--------|-------------|
| `gvcf.sh` | Walkthrough: generate a gVCF from a BAM using GATK HaplotypeCaller |
| `validate.sh` | Validate a gVCF with GATK ValidateVariants |
| `annotate_clinvar.py` | Match gVCF variants against ClinVar; outputs TSV + markdown report |
| `roh.sh` | Runs-of-homozygosity analysis using bcftools + PLINK |
| `mito_haplogroup.py` | Candidate mitochondrial haplogroup screen (limited local PhyloTree subset) |
| `mito_heteroplasmy.py` | Mitochondrial heteroplasmy landscape: plots, age hotspots, MITOMAP screen |
| `pca/01_download.sh` | Download 1000 Genomes Phase 3 reference data (~12.7 GB) |
| `pca/02_pca.sh` | Ancestry PCA: project your genome onto 2504 reference samples |
| `pca/03_plot.py` | Plot PCA results (PC1/PC2, PC1/PC3) |

## Requirements

- `bcftools`, `samtools` (htslib)
- GATK 4.6+ (`gvcf.sh` downloads it)
- `plink2` and PLINK 1.x (`roh.sh`, `pca/`)
- Python 3 with packages in `requirements.txt`

## License

0BSD. See [LICENSE](LICENSE).
