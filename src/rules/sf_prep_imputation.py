"""Make white ppl vcf files for imputation server."""

rule restrict_white:
    input:
        b = DATA + 'interim/bfiles_filter_samples/{group}.fam',
        k = DATA + 'interim/mds_cut/{group}.keep_samples'
    output:
        DATA + 'interim/bfiles_eur/{group}.fam',
    singularity:
        PLINK
    log:
        LOG + 'prs/restrict.eur.{group}'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_samples/{wildcards.group} --keep {input.k} "
        "--make-bed --out {DATA}interim/bfiles_eur/{wildcards.group} &> {log}"

rule vcf_white:
    input:
        DATA + 'interim/bfiles_eur/{group}.fam',
    output:
        DATA + 'interim/bfiles_eur_vcf_chr{chr}/{group}.vcf'
    singularity:
        PLINK
    log:
        LOG + 'prs/vcf.{chr}.{group}.eur'
    shell:
        "plink --bfile {DATA}interim/bfiles_eur/{wildcards.group} --recode vcf --chr {wildcards.chr} "
        "--out {DATA}interim/bfiles_eur_vcf_chr{wildcards.chr}/{wildcards.group} &> {log}"

rule hack_tpop_fam:
    input:
        DATA + 'interim/bfiles_filter_samples/{group}.fam',
    output:
        DATA + 'interim/bfiles_tpop/{group}.fam',
    shell:
        'cp {input} {output}'

rule vcf_tpop:
    input:
        DATA + 'interim/bfiles_filter_samples/{group}.fam',
    output:
        DATA + 'interim/bfiles_tpop_vcf_chr{chr}/{group}.vcf'
    singularity:
        PLINK
    log:
        LOG + 'prs/vcf.{chr}.{group}.tpop'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_samples/{wildcards.group} "
        "--recode vcf --chr {wildcards.chr} "
        "--out {DATA}interim/bfiles_tpop_vcf_chr{wildcards.chr}/{wildcards.group} &> {log}"

rule cp_vcf:
    input:
        DATA + 'interim/bfiles_{pop}_vcf_chr{chr}/{group}.vcf'
    output:
        DATA + 'interim/{pop}_vcf/{group}.{chr}.vcf'
    shell:
        'cp {input} {output}'

rule zip_vcf:
    input:
        DATA + 'interim/{pop}_vcf/{group}.{chr}.vcf'
    output:
        DATA + "interim/{pop}_vcf/{group}.{chr}.vcf.gz"
    wrapper:
        "0.27.1/bio/vcf/compress"
