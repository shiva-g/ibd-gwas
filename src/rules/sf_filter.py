rule filter_snps:
    """rm snps missing in more than 5% of samples: --geno
       rm non-polymorphic: --min-ac
    """
    input:
        expand(DATA + 'interim/bfiles/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/bfiles_filter_snps_1/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.filter_snps'
    shell:
        "plink --bfile {DATA}interim/bfiles/{wildcards.group} --min-ac 1 --geno 0.05 "
        "--make-bed --out {DATA}interim/bfiles_filter_snps_1/{wildcards.group} &> {log}"

rule list_dup_pos:
    input:  i = DATA + 'interim/bfiles_filter_snps_1/{group}.bim'
    output: o = DATA + 'interim/bfiles_filter_snps_1_dupPos/{group}.duppos'
    run:
        names = ['chrom', 'id', 'blank', 'pos', 'allele1', 'allele2']
        df = pd.read_csv(input.i, header=None, sep='\t', names = names)
        size_df = df[['pos']].groupby('pos').size().reset_index().rename(columns={0:'s'})
        m = pd.merge(size_df[size_df.s>1], df, how='left', on='pos')
        m[['id']].to_csv(output.o, index=False, header=None)

rule list_dup_snps:
    input:
        expand(DATA + 'interim/bfiles_filter_snps_1/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/bfiles_filter_snps_dups/{group}.dupvar'
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.dup_snps'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_snps_1/{wildcards.group} --list-duplicate-vars ids-only "
        "--out {DATA}interim/bfiles_filter_snps_dups/{wildcards.group} &> {log}"

rule rm_dup_snps:
    input:
        b=expand(DATA + 'interim/bfiles_filter_snps_1/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') ),
        d=DATA + 'interim/bfiles_filter_snps_dups/{group}.dupvar',
        dd = DATA + 'interim/bfiles_filter_snps_1_dupPos/{group}.duppos'
    output:
        expand(DATA + 'interim/bfiles_filter_snps/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.rm_dup_snps'
    shell:
        'plink --bfile {DATA}interim/bfiles_filter_snps_1/{wildcards.group} '
        '--exclude <(cat {input.dd} <(sed "s/ /\\n/g" {input.d}) | sort -u) '
        "--make-bed --out {DATA}interim/bfiles_filter_snps/{wildcards.group} &> {log}"

rule filter_samples:
    """rm samples missing in more than 5% of snps: --mind
    """
    input:
        expand(DATA + 'interim/bfiles_filter_snps/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/bfiles_filter_samples/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.filter_samples'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_snps/{wildcards.group} --mind 0.05 --make-bed "
        "--out {DATA}interim/bfiles_filter_samples/{wildcards.group} &> {log}"

rule indep_snps:
    input:
        expand(DATA + 'interim/bfiles_filter_samples/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/bfiles_indep/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.indep_snps'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_samples/{wildcards.group} --indep 50 5 2 "
        "--make-bed --out {DATA}interim/bfiles_indep/{wildcards.group} &> {log}"

rule missing:
    input:
        expand(DATA + 'interim/bfiles_filter_snps/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/missing/{{group}}.{miss}', miss=('imiss', 'lmiss') )
    singularity:
        PLINK
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_snps/{wildcards.group} --missing --out {DATA}interim/missing/{wildcards.group}"
