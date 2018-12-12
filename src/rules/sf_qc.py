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
        expand(DATA + 'interim/bfiles_filter_snps_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/missing_test/{{group}}.{miss}', miss=('imiss', 'lmiss') )
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.missing'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_snps_nox/{wildcards.group} "
        "--missing --test-missing --out {DATA}interim/missing_test/{wildcards.group} &> {log}"

rule check_freq_before_imputation:
    input:
        f = DATA + 'interim/bfiles_filter_samples_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_filter_samples_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/qc_freq_before_impute/{group}.frq'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.freq'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} "
        "--freq --out $(dirname {output})/{wildcards.group} &> {log}"

rule check_freq_after_imputation:
    input:
        f = DATA + 'interim/bfiles_filter_samples_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_filter_samples_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/qc_freq_after_impute/{group}.frq'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.freq'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} "
        "--freq --out $(dirname {output})/{wildcards.group} &> {log}"

rule summarize_freq:
    input:
        i = DATA + 'interim/qc_freq_{imputeStatus}_impute/{group}.frq'
    output:
        o = DATA + 'interim/qc_freq_{imputeStatus}_impute/{group}.counts'
    run:
        df = pd.read_csv(input.i, delim_whitespace=True)
        def assign_class(row):
            if row['MAF']<0.01:
                return 'Below 1%'
            elif row['MAF']>=0.01 and row['MAF']<0.05:
                return '1% to 5%'
            else:
                return 'Greater 5%'

        df.loc[:, 'group'] = df.apply(assign_class, axis=1)
        df.groupby('group').size().reset_index().to_csv(output.o, index=False, sep='\t')

rule check_hwe:
    input:
        expand(DATA + 'interim/bfiles_filter_samples_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/qc_hwe/{group}.hwe'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.hwe'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_samples_nox/{wildcards.group} "
        "--hardy --out {DATA}interim/qc_hwe/{wildcards.group} &> {log}"

rule summarize_hwe:
    input:
        i = DATA + 'interim/qc_hwe/{group}.hwe'
    output:
        o = DATA + 'interim/qc_hwe/{group}.counts'
    run:
        df = pd.read_csv(input.i, delim_whitespace=True)
        crit = df.apply(lambda row: row['TEST'] in ('AFF', 'UNAFF'), axis=1)
        df = df[crit]
        crit01 = df.apply(lambda row: row['P']<0.01, axis=1)
        crit001 = df.apply(lambda row: row['P']<0.001, axis=1)
        cols = {0:'count'}
        g01 = df[crit01].groupby('TEST').size().reset_index().rename(columns=cols)
        g01.loc[:, 'pCut'] = 0.01
        g001 = df[crit001].groupby('TEST').size().reset_index().rename(columns=cols)
        g001.loc[:, 'pCut'] = 0.001
        pd.concat([g01, g001]).to_csv(output.o, index=False, sep='\t')

rule all_qc:
    input: expand(DATA + 'interim/missing_test/3groups.{miss}', miss=('imiss', 'lmiss') ),
           DATA + 'interim/sex_check/3groups.sexcheck',
           DATA + 'interim/qc_hwe/3groups.counts',
           DATA + 'interim/qc_freq/3groups.counts'
