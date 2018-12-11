"""Main snakefile for polygeic risk score"""

include: "const.py"
include: "sf_manifest.py"
include: "sf_filter.py"
include: "sf_prep.py"
include: "sf_qc.py"
#include: "https://raw.githubusercontent.com/samesense/snakemake-liftover-workflow/master/Snakefile.py"
include: "sf_hapmap.py"
include: "sf_mds.py"
include: "sf_prs.py"

rule before_imputation:
    input:
        expand(DATA + 'interim/eur_vcf/3groups.{chr}.vcf.gz', chr=range(1,23)),
        PLOTS + 'eur_mds.png',
        PLOTS + 'hapmap_mds.png',
        expand(DATA + 'interim/missing_test/3groups.{miss}', miss=('imiss', 'lmiss') ),
        DATA + 'interim/sex_check/3groups.sexcheck',
        DATA + 'interim/qc_hwe/3groups.counts',
        DATA + 'interim/qc_freq/3groups.counts'

rule after_imputation:
    input:
        DATA + 'interim/prsice/eur.summary'

