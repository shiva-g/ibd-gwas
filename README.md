ibd-gwas
==============================

IBD vs healthy
* Manifest: /mnt/isilon/microbiome/perry/ibd-gwas/data/processed/MANIFEST.csv
* [Methods](writeup/methods.md)
* [Log](writeup/log.md)
* run: 
`conda create env create -f conda_reqs.yml
# before imputation
snakemake -s Snakefile.py --use-singularity --singularity-args "-B /mnt/isilon/:/mnt/isilon" --use-conda -j22 before_imputation
# impute with https://imputationserver.sph.umich.edu/index.html#!
# download imputed files to data/interim/imputed/
# run after imputation
snakemake -s Snakefile.py --use-singularity --singularity-args "-B /mnt/isilon/:/mnt/isilon" --use-conda -j22 after_imputation
`
