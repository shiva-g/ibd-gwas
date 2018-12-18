ibd-gwas
==============================

IBD vs healthy
* Manifest: /mnt/isilon/microbiome/perry/ibd-gwas/data/processed/MANIFEST.csv
* [Results](writeup/results.md)
* [Methods](writeup/methods.md)
* [Log](writeup/log.md)
* [Meetings](writeup/meetings.md)

### Run: 
`conda create env create -f envs/project.yml`
#### before imputation
`snakemake -s Snakefile.py --use-singularity --singularity-args "-B /mnt/isilon/:/mnt/isilon" --use-conda -j22 before_imputation`
#### impute with https://imputationserver.sph.umich.edu/index.html w/ data/iterim/eur_vcf/*vcf.gz and download imputed files to data/interim/imputed/
#### run after imputation
`snakemake -s Snakefile.py --use-singularity --singularity-args "-B /mnt/isilon/:/mnt/isilon" --use-conda -j22 after_imputation
`
