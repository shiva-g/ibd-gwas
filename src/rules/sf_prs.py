"""Polygenic risk score"""

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
        LOG + 'prs/vcf.{chr}.{group}'
    shell:
        "plink --bfile {DATA}interim/bfiles_eur/{wildcards.group} --recode vcf --chr {wildcards.chr} "
        "--out {DATA}interim/bfiles_eur_vcf_chr{wildcards.chr}/{wildcards.group} &> {log}"

rule cp_vcf:
    input:
        DATA + 'interim/bfiles_eur_vcf_chr{chr}/{group}.vcf'
    output:
        DATA + 'interim/eur_vcf/{group}.{chr}.vcf'
    shell:
        'cp {input} {output}'

rule zip_vcf:
    input:
        DATA + 'interim/eur_vcf/{group}.{chr}.vcf'
    output:
        DATA + "interim/eur_vcf/{group}.{chr}.vcf.gz"
    wrapper:
        "0.27.1/bio/vcf/compress"

# https://www.nature.com/articles/ng.3359#supplementary-information
rule prep_gwas:
    input:
        i = DATA + 'raw/ibd_gwas/ng.3359-S4.xlsx'
    output:
        o = DATA + 'interim/ibd_gwas.assoc'
    run:
        df = pd.read_excel(input.i, sheet_name='Heterogeneity of effect', skiprows=7)
        cols = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'EUR_OR', 'EUR_PVAL', 'EUR_SE']
        df[cols].to_csv(output.o, index=False, sep=' ')

rule for_imputation:
    input: expand(DATA + 'interim/eur_vcf/3groups.{chr}.vcf.gz', chr=range(1,23))

rule unzip_imputed:
    input:
        DATA + 'interim/imputed/chr_{c}.zip'
    output:
        DATA + 'interim/imputed/chr{c}.dose.vcf.gz'
    shell:
        'unzip -o -P VUekgAD02lfhM5 -d {DATA}/interim/imputed/ {input}'

rule vcf_to_bfiles:
    input:
        DATA + 'interim/imputed/chr{chr}.dose.vcf.gz'
    output:
        expand(DATA + 'interim/bfiles_imputed/chr{{chr}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prs/{chr}.vcf'
    shell:
        "plink --vcf {input} --make-bed --const-fid 0 "
        "--out {DATA}interim/bfiles_imputed/chr{wildcards.chr} &> {log}"

rule merge_imputed_bfiles_ls:
    input:
        expand(DATA + 'interim/bfiles_imputed/chr{chr}.fam', chr=range(1, 23))
    output:
        o = DATA + 'interim/bfiles_imputed.eur.ls'
    run:
        with open(output.o, 'w') as fout:
            for afile in list(input):
                print(afile.split('.')[0], file=fout)

rule merge_imputed_bfiles:
    input:
        DATA + 'interim/bfiles_imputed.eur.ls'
    output:
        DATA + 'processed/bfiles_imputed/eur.fam'
    singularity:
        PLINK
    log:
        LOG + 'prs/merge'
    shell:
        "plink --merge-list {input} --make-bed "
        "--out {DATA}interim/bfiles_imputed_tmp/eur &> {log}.1 || "
        "plink --merge-list {input} --biallelic-only strict "
        "--exclude {DATA}interim/bfiles_imputed_tmp/eur-merge.missnp --make-bed "
        "--out {DATA}processed/bfiles_imputed/eur &> {log} "

rule prsice:
    input:
        b = DATA + 'processed/bfiles_imputed/eur.fam',
        a = DATA + 'interim/ibd_gwas.assoc'
    output:
        DATA + 'interim/prsice/eur'
    singularity:
        PRSICE
    threads: 8
    shell:
        'Rscript /usr/local/bin/PRSice.R --dir {DATA}/interim/prsice/ '
        '--snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat EUR_OR --se EUR_SE --pvalue EUR_PVAL '
        '--prsice /usr/local/bin/PRSice_linux '
        '--base {input.a} '
        '--target {DATA}processed/bfiles_imputed/eur '
        '--thread {threads} --stat OR --binary-target T '
        '--out {DATA}interim/prsice/eur'

rule unzip_imputation:
    input: expand(DATA + 'interim/imputed/chr{chr}.dose.vcf.gz', chr=range(1,23))
