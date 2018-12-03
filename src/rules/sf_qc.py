rule check_gender_y:
    input:
        expand(DATA + 'interim/bfiles_indep/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/sex_check_y/{group}.sexcheck'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.sex_check_y'
    shell:
        "plink --bfile {DATA}interim/bfiles_indep/{wildcards.group} --check-sex y-only "
        "--out {DATA}interim/sex_check_y/{wildcards.group} &> {log}"

rule check_gender_x:
    input:
        expand(DATA + 'interim/bfiles_indep/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/sex_check/{group}.sexcheck'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.sex_check'
    shell:
        "plink --bfile {DATA}interim/bfiles_indep/{wildcards.group} --check-sex "
        "--out {DATA}interim/sex_check/{wildcards.group} &> {log}"

rule missing:
    input:
        expand(DATA + 'interim/bfiles_filter_snps/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/missing/{{group}}.{miss}', miss=('imiss', 'lmiss') )
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.missing'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_snps/{wildcards.group} --missing --out {DATA}interim/missing/{wildcards.group} &> {log}"

rule all_qc:
    input: expand(DATA + 'interim/missing/3groups.{miss}', miss=('imiss', 'lmiss') ),
           DATA + 'interim/sex_check/3groups.sexcheck',
#           DATA + 'interim/sex_check_y/3groups.sexcheck'
