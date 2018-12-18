"""Annotate imputed vcf file"""

rule vcf_imputed:
    input:
        DATA + 'processed/bfiles_imputed/{group}.fam',
    output:
        DATA + 'interim/bfiles_clean_vcf/{group}.vcf'
    singularity:
        PLINK
    log:
        LOG + 'ann/vcf.{group}'
    shell:
        "plink --bfile {DATA}processed/bfiles_imputed/{wildcards.group} --recode vcf "
        "--out {DATA}interim/bfiles_clean_vcf/{wildcards.group} &> {log}"

rule snpeff:
    input:
        DATA + "interim/bfiles_clean_vcf/{group}.vcf"
    output:
        vcf=DATA + "interim/variants/snpeff/{group}.vcf",    # the main output file, required
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="GRCh37.75", # reference name (from `snpeff databases`)
        extra="-Xmx32g -Xms16g"  # optional parameters (e.g., max memory 4g)
    singularity:
        "docker://quay.research.chop.edu/evansj/snpeff-docker:grch37"
    conda:
        ENVS + "project.yml"
    threads: 30
    #
    #shell:
    #    "snpEff -noStats -t {threads} -Xmx32g -Xms16g {params.reference} {input} > {output} 2> {log}"
    #
    wrapper:
        "file:///mnt/isilon/cbmi/variome/perry/projects/ext/snakemake-wrappers/bio/snpeff/wrapper.py"
