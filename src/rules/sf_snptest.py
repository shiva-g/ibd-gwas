"""Run snptest"""

# https://swvanderlaan.github.io/post/converting-m3vcf-to-anything-else/
rule gen:
    input:
        DATA + 'interim/imputed_r2_limit_vcf/chr{c}.vcf.gz'
    output:
        o = DATA + 'interim/eur_gen/chr{c}.tmp.gen'
    singularity:
        PLINK2
    threads:
        32
    shell:
        "plink2 --vcf {input} dosage=GP --export oxford --const-fid 0 "
        "--threads {threads} --out {DATA}interim/eur_gen/chr{wildcards.c}.tmp"

rule update_gen_sample:
    """Mk fam for imputed files b/c imputation removed data.
       Some data will be missing for samples removed from study.
    """
    input:
        gen = DATA + 'interim/eur_gen/chr{c}.tmp.gen',
        data_fam = DATA + 'interim/bfiles_eur/3groups.fam'
    output:
        o = DATA + 'interim/eur_gen/chr{c}.sample',
        gen = DATA + 'interim/eur_gen/chr{c}.gen'
    run:
        sample_file = input.gen.replace('.gen', '.sample')
        use_cols = ['ID_1', 'ID_2', 'missing', 'sex', 'pheno']
        use_cols_short = ['ID_1', 'sex', 'pheno']
        cols= ['fid', 'iid', 'f', 'm', 'sex', 'pheno']
        int_cols= ['fid', 'f', 'm', 'sex', 'pheno']
        dtype={'fid':int, 'f':int, 'm':int, 'sex':int, 'pheno':int}
        df_sample = pd.read_csv(sample_file, delim_whitespace=True)[['ID_1', 'ID_2', 'missing']]
        df_sample.loc[:, 'ID_1'] = df_sample['ID_2']
        df_sample['ID_2'] = 0
        df_dat = pd.read_csv(input.data_fam, sep=' ', header=None, names=cols, dtype=dtype)
        df_dat.loc[:, 'ID_1'] = df_dat.apply(lambda row: '0_' + row['iid'], axis=1)
        # failes for NAs
        df_dat.loc[:, 'pheno'] = df_dat.apply(lambda row: '0' if row['pheno']==1 else '1', axis=1)
        bad_samples = set(df_sample['ID_1']) - set(df_dat['ID_1'])
        bad_dat = []
        for bad_sample in bad_samples:
            if bad_sample == '0':
                bad_dat.append( [bad_sample, 'D', 'B',])
            else:
                bad_dat.append( [bad_sample, '0', 'NA',])

        bad_df = pd.DataFrame(bad_dat, columns=use_cols_short)
        df = pd.merge(df_sample, pd.concat([df_dat[['ID_1', 'sex', 'pheno']], bad_df]), how='left', on='ID_1')
        #df.loc[:, 'ID_1'] = df.apply(lambda row: row['ID_1'].replace('_','a'), axis=1)
        df[use_cols].to_csv(output.o, sep=' ', index=False)
        shell('cp {input.gen} {output.gen}')

# association
rule snptest:
    input:
        gen = DATA + 'interim/eur_gen/chr{c}.gen',
        sample  = DATA + 'interim/eur_gen/chr{c}.sample',
    output:
        DATA + 'interim/eur_snptest/chr{c}.out'
    log:
        LOG + 'snptest/eur.chr{c}'
    singularity:
        SNPTEST
    shell:
        """
        snptest \
        -data {input.gen} {input.sample} \
        -o {output} \
        -genotype_field GP \
        -frequentist 1 \
        -method score \
        -pheno pheno \
        -hwe
        """

rule gens:
    input:
        expand(DATA + 'interim/eur_snptest/chr{c}.out', c=range(22,23))
