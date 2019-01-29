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

rule check_hwe:
    input:
        f = DATA + 'interim/bfiles_filter_samples_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_filter_samples_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/qc_hwe/{group}.hwe'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.hwe'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} "
        "--hardy --out $(dirname {output})/{wildcards.group} &> {log}"

rule check_het:
    input:
        f = DATA + 'interim/bfiles_indep_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_indep_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/qc_het/{group}.het'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.het'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} "
        "--het --out $(dirname {output})/{wildcards.group} &> {log}"

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
    input: expand(DATA + 'interim/missing_test/{group}.{miss}', group=('3groups', 'gsa'), miss=('imiss', 'lmiss') ),
           # expand(DATA + 'interim/sex_check/{group}.sexcheck', group=('3groups', 'gsa')),
           expand(DATA + 'interim/qc_hwe/{group}.counts', group=('3groups', 'gsa')),
           # expand(DATA + 'interim/qc_freq/{group}.counts', group=('3groups', 'gsa')),
           expand(DATA + 'interim/plink_genome/{group}.genome.tab', group=('3groups', 'gsa')),
           expand(DATA + 'interim/qc_het/{group}.het', group=('3groups', 'gsa')),
