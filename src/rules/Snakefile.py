"""Main snakefile for polygeic risk score"""

include: "const.py"
include: "sf_manifest.py"
include: "sf_mafs.py"
include: "sf_filter.py"
include: "sf_prep.py"
include: "sf_prep_veo.py"
include: "sf_qc.py"
include: "sf_liftover.py"
#include: "https://raw.githubusercontent.com/samesense/snakemake-liftover-workflow/master/Snakefile.py"
include: "sf_hapmap.py"
include: "sf_mds.py"
include: "sf_prep_imputation.py"
include: "sf_clean_imputed.py"
include: "sf_prs.py"
include: "sf_assoc.py"
include: "sf_snptest.py"
include: "sf_ann.py"
include: "sf_snp_splits.py"

rule before_imputation:
    input:
        expand(DATA + 'interim/{pop}_vcf/3groups.{chr}.vcf.gz', pop=('tpop', 'eur'), chr=range(1,23)),
        PLOTS + 'eur_mds.png',
        PLOTS + 'hapmap_mds.png',
        PLOTS + 'hapmap_mds_nogroup.png',
        expand(DATA + 'interim/missing_test/3groups.{miss}', miss=('imiss', 'lmiss') ),
        DATA + 'interim/sex_check/3groups.sexcheck',
        DATA + 'interim/qc_hwe/3groups.counts',
        DATA + 'interim/qc_het/3groups.het'

rule after_imputation:
    input:
        expand(PLOTS + '{group}.{pop}.{go}.prs.roc.png', group=G, pop=('eur', 'tpop'), go=GOS),
        expand(PLOTS + '{group}.{pop}.{go}.prs.density.png', group=G, pop=('eur', 'tpop'), go=GOS),
        #expand(PLOTS + '{group}.{pop}.{go}.prs.quantiles.png', group=G, pop=('eur', 'tpop'), go=GOS),
        # # expand(DATA + 'interim/plink_assoc_fmt/{group}/eur.assoc', group=G),
        # # expand(PLOTS + 'manhattan.{group}.png', group=G),
        expand(PWD + 'writeup/tables/prs.{pop}.md', pop=('eur', 'tpop')),
        # # PWD + 'writeup/tables/maf.md',
        # # expand(DATA + 'interim/prsice/snp_overlap/all.{pop}.init', pop=('eur', 'tpop')),
        #expand(PWD + "writeup/tables/{age}.{group}.{pop}.assoc.csv", age=('ped', 'adult'), pop=('tpop', 'eur'), group=G),
        # DATA + 'interim/snp_groups/tpop'
