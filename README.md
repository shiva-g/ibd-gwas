ibd-gwas
==============================

IBD vs healthy
* Manifest: /mnt/isilon/microbiome/perry/ibd-gwas/data/processed/MANIFEST.csv
* [Methods](writeup/methods.md)
* [Log](writeup/log.md)
* run: `snakemake -s Snakefile.py --use-singularity --singularity-args "-B /mnt/isilon/:/mnt/isilon" --use-conda -j22 prs`
