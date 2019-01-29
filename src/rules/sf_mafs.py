"""Rm chr x and calc mafs."""

rule drop_x_filter_samples:
    """rm x for hwe test
    """
    input:
        f = DATA + 'interim/bfiles_filter_samples/{group}.fam',
        b = expand(DATA + 'interim/bfiles_filter_samples/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        f = DATA + 'interim/bfiles_filter_samples_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_filter_samples_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.filter_samples_nox'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} --chr 1-22 --make-bed "
        "--out $(dirname {output.f})/{wildcards.group} &> {log}"

rule drop_x_imputed_prep:
    """rm x for hwe test
    """
    input:
        f = DATA + 'interim/bfiles_{pop}/{group}.fam',
        b = expand(DATA + 'interim/bfiles_{{pop}}/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        f = DATA + 'interim/bfiles_{pop,tpop|eur}_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_{{pop,tpop|eur}}_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/{group}.{pop}_nox'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} --chr 1-22 --make-bed "
        "--out $(dirname {output.f})/{wildcards.group} &> {log}"

# rule drop_x_imputed:
#     input:
#         f = DATA + 'interim/bfiles_imputed_combined/{group}.fam',
#         b = expand(DATA + 'interim/bfiles_imputed_combined/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     output:
#         f = DATA + 'interim/bfiles_imputed_{pop}_nox/{group}.fam',
#         b = expand(DATA + 'interim/bfiles_imputed_{{pop}}_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     singularity:
#         PLINK
#     log:
#         LOG + 'prep/{group}.{pop}_nox'
#     shell:
#         "plink --bfile $(dirname {input.f})/{wildcards.group} --chr 1-22 --make-bed "
#         "--out $(dirname {output.f})/{wildcards.group} &> {log}"

# rule drop_x_indep:
#     """rm x for freq test
#     """
#     input:
#         f = DATA + 'interim/bfiles_indep/{group}.fam',
#         b = expand(DATA + 'interim/bfiles_indep/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     output:
#         f = DATA + 'interim/bfiles_indep_nox/{group}.fam',
#         b = expand(DATA + 'interim/bfiles_indep_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     singularity:
#         PLINK
#     log:
#         LOG + 'prep/{group}.indep_nox'
#     shell:
#         "plink --bfile $(dirname {input.f})/{wildcards.group} --chr 1-22 --make-bed "
#         "--out $(dirname {output.f})/{wildcards.group} &> {log}"

rule check_freq_notfor_imputation:
    input:
        f = DATA + 'interim/bfiles_filter_samples_nox/{group}.fam',
        b = expand(DATA + 'interim/bfiles_filter_samples_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        DATA + 'interim/qc_freq_notfor_impute/{group}.frq'
    singularity:
        PLINK
    log:
        LOG + 'qc/{group}.freq'
    shell:
        "plink --bfile $(dirname {input.f})/{wildcards.group} "
        "--freq --out $(dirname {output})/{wildcards.group} &> {log}"

#rule check_freq_usedfor_imputation:
#     input:
#         f = DATA + 'interim/bfiles_{pop}_nox/{group}.fam',
#         b = expand(DATA + 'interim/bfiles_{{pop}}_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     output:
#         DATA + 'interim/qc_freq_usedfor_impute/{group}.{pop}.frq'
#     singularity:
#         PLINK
#     log:
#         LOG + 'qc/{group}.{pop}.used_for.freq'
#     shell:
#         "plink --bfile $(dirname {input.f})/{wildcards.group} "
#         "--freq --out $(dirname {output})/{wildcards.group} &> {log}"

# rule check_freq_after_imputation:
#     input:
#         f = DATA + 'interim/bfiles_imputed_{pop}_nox/{group}.fam',
#         b = expand(DATA + 'interim/bfiles_imputed_{{pop}}_nox/{{group}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     output:
#         DATA + 'interim/qc_freq_after_impute/{group}.{pop}.frq'
#     singularity:
#         PLINK
#     log:
#         LOG + 'qc/{group}.after_impute.freq'
#     shell:
#         "plink --bfile $(dirname {input.f})/{wildcards.group} "
#         "--freq --out $(dirname {output})/{wildcards.group} &> {log}"

#rule fake_gsa_after_impute:
 #    input:
 #        DATA + 'interim/qc_freq_usedfor_impute/gsa.frq'
 #    output:
 #        DATA + 'interim/qc_freq_after_impute/gsa.tpop.frq'
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

# rule combine_mafs:
#     input:
#         qc=DATA + 'interim/qc_freq_notfor_impute/{group}.counts',
#         eur=DATA + 'interim/qc_freq_usedfor_impute/{group}.{pop}.counts',
#         imp=DATA + 'interim/qc_freq_after_impute/{pop}.counts',
#     output:
#         o = PWD + 'writeup/tables/{group}.{pop}.maf.md'
#     run:
#         def read_df(afile, label):
#             df = pd.read_csv(afile, sep='\t').rename(columns={'0':'snp_count'})
#             df['dataset'] = label
#             return df
#         with open(output.o, 'w') as fout:
#             df = pd.concat([read_df(input.qc, 'afterQC'), read_df(input.eur, 'beforeImpute'), read_df(input.imp, 'afterImpute')])
#             print(tabulate.tabulate(df.values, df.columns, tablefmt="pipe"), file=fout)

rule combine_mafs:
    input:
        qc=DATA + 'interim/qc_freq_notfor_impute/{group}.counts',
    output:
        o = PWD + 'writeup/tables/{group}.maf.md'
    run:
        def read_df(afile, label):
            df = pd.read_csv(afile, sep='\t').rename(columns={'0':'snp_count'})
            df['dataset'] = label
            return df
        with open(output.o, 'w') as fout:
            df = read_df(input.qc, 'afterQC')
            print(tabulate.tabulate(df.values, df.columns, tablefmt="pipe"), file=fout)

rule maf_tables:
    input:
        PWD + 'writeup/tables/gsa.maf.md'
