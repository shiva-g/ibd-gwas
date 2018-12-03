rule check_gender_x:
    input:
        expand(DATA + 'interim/bfiles/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/sex_check/{group}.x'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.sex_check'
    shell:
        "plink --bfile {DATA}interim/bfiles/{wildcards.group} --check-sex "
        "--out {DATA}interim/sex_check/{wildcards.group} &> {log}"

rule missing:
    input:
        expand(DATA + 'interim/bfiles_filter_snps/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/missing/{{group}}.{miss}', miss=('imiss', 'lmiss') )
    singularity:
        PLINK
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_snps/{wildcards.group} --missing --out {DATA}interim/missing/{wildcards.group}"

rule all_qc:
    input: expand(DATA + 'interim/missing/3groups.{miss}', miss=('imiss', 'lmiss') )
