"""Main snakefile for polygeic risk score"""

include: "const.py"
include: "sf_manifest.py"
include: "sf_mafs.py"
include: "sf_filter.py"
include: "sf_prep.py"
include: "sf_qc.py"
include: "sf_liftover.py"
#include: "https://raw.githubusercontent.com/samesense/snakemake-liftover-workflow/master/Snakefile.py"
include: "sf_hapmap.py"
include: "sf_mds.py"
include: "sf_prep_imputation.py"
include: "sf_clean_imputed.py"
include: "sf_prs.py"
include: "sf_assoc.py"
include: "sf_ann.py"

rule before_imputation:
    input:
        expand(DATA + 'interim/eur_vcf/3groups.{chr}.vcf.gz', chr=range(1,23)),
        PLOTS + 'eur_mds.png',
        PLOTS + 'hapmap_mds.png',
        PLOTS + 'hapmap_mds_nogroup.png',
        expand(DATA + 'interim/missing_test/3groups.{miss}', miss=('imiss', 'lmiss') ),
        DATA + 'interim/sex_check/3groups.sexcheck',
        DATA + 'interim/qc_hwe/3groups.counts',
        DATA + 'interim/qc_het/3groups.het'

rule after_imputation:
    input:
        expand(PLOTS + '{group}.eur.prs.roc.png', group=G),
        expand(PLOTS + '{group}.eur.prs.density.png', group=G),
        expand(DATA + 'interim/plink_assoc_fmt/{group}/eur.assoc', group=G),
        expand(PLOTS + 'manhattan.{group}.png', group=G),
        PWD + 'writeup/tables/prs.eur.md',
        PWD + 'writeup/tables/maf.md',
        DATA + "interim/variants/snpeff/eur.vcf",
        DATA + 'interim/prsice/snp_overlap/all.eur.init'
