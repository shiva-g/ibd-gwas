"""Run snptest"""

# https://swvanderlaan.github.io/post/converting-m3vcf-to-anything-else/
rule gen:
    input:
        DATA + 'interim/imputed_r2_limit_vcf_{pop}/chr{c}.vcf.gz'
    output:
        o = DATA + 'interim/{pop}_gen/chr{c}.tmp.gen'
    singularity:
        PLINK2
    threads:
        32
    shell:
        "plink2 --vcf {input} dosage=GP --export oxford --const-fid 0 "
        "--threads {threads} --out {DATA}interim/{wildcards.pop}_gen/chr{wildcards.c}.tmp"

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
        df[use_cols].to_csv(output.o, sep=' ', index=False)
        shell('cp {input.gen} {output.gen}')

rule update_gen_sample_tpop:
    """Mk fam for imputed files b/c imputation removed data.
       Some data will be missing for samples removed from study.
    """
    input:
        gen = DATA + 'interim/tpop_gen/chr{c}.tmp.gen',
        data_fam = DATA + 'interim/bfiles_tpop/3groups.fam',
        pcs = DATA + 'interim/mds_dat/ibd_hapmap.dat'
    output:
        o = DATA + 'interim/tpop_gen/chr{c}.sample',
        gen = DATA + 'interim/tpop_gen/chr{c}.gen'
    run:
        pcs = pd.read_csv(input.pcs, sep='\t')
        sample_file = input.gen.replace('.gen', '.sample')
        use_cols = ['ID_1', 'ID_2', 'missing', 'sex', 'pheno', 'C1', 'C2']
        use_cols_short = ['ID_1', 'sex', 'pheno', 'C1', 'C2']
        cols= ['fid', 'iid', 'f', 'm', 'sex', 'pheno']
        int_cols= ['fid', 'f', 'm', 'sex', 'pheno']
        dtype={'fid':int, 'f':int, 'm':int, 'sex':int, 'pheno':int}
        df_sample = pd.read_csv(sample_file, delim_whitespace=True)[['ID_1', 'ID_2', 'missing']]
        df_sample.loc[:, 'ID_1'] = df_sample['ID_2']
        df_sample['ID_2'] = 0
        df_dat = pd.read_csv(input.data_fam, sep=' ', header=None, names=cols, dtype=dtype)
        df_dat.loc[:, 'ID_1'] = df_dat.apply(lambda row: '0_' + row['iid'], axis=1)
        pcs.loc[:, 'ID_1'] = pcs.apply(lambda row: '0_' + row['IID'], axis=1)
        pcs = pcs[['ID_1', 'C1', 'C2']]
        df_dat = pd.merge(df_dat, pcs, on='ID_1', how='left')
        # fails for NAs
        df_dat.loc[:, 'pheno'] = df_dat.apply(lambda row: '0' if row['pheno']==1 else '1', axis=1)
        bad_samples = set(df_sample['ID_1']) - set(df_dat['ID_1'])
        bad_dat = []
        for bad_sample in bad_samples:
            if bad_sample == '0':
                bad_dat.append( [bad_sample, 'D', 'B', 'C', 'C'])
            else:
                bad_dat.append( [bad_sample, '0', 'NA', 'NA', 'NA'])

        bad_df = pd.DataFrame(bad_dat, columns=use_cols_short)
        df = pd.merge(df_sample, pd.concat([df_dat[['ID_1', 'sex', 'pheno', 'C1', 'C2']], bad_df]), how='left', on='ID_1')
        df[use_cols].to_csv(output.o, sep=' ', index=False)
        shell('cp {input.gen} {output.gen}')

rule tmp_p:
    input:
        DATA + 'interim/tpop_snptest/chr22.out',

# association eur
rule snptest_eur:
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
        -hwe &> {log}
        """

# association correct for pop structure
rule snptest_tpop:
    input:
        gen = DATA + 'interim/tpop_gen/chr{c}.gen',
        sample  = DATA + 'interim/tpop_gen/chr{c}.sample',
    output:
        DATA + 'interim/tpop_snptest/chr{c}.out'
    log:
        LOG + 'snptest/tpop.chr{c}'
    singularity:
        SNPTEST
    shell:
        "snptest -data {input.gen} {input.sample} "
        "-o {output} -genotype_field GP -frequentist 1 "
        "-method score -pheno pheno -cov_names C1 C2 -hwe &> {log}"

rule snptest_collapse:
    input:
        expand(DATA + 'interim/{{pop}}_snptest/chr{c}.out', c=range(1, 23))
    output:
        o = DATA + 'interim/{pop}_snptest_final/snptest.out'
    run:
        def read_csv(afile):
            chrom = afile.split('/')[-1].split('.')[0][3:]
            df = pd.read_csv(afile, comment='#', sep=' ').drop(['chromosome',], axis=1)
            df = df.assign(chromosome=chrom)
            return df

        dfs = [read_csv(_) for _ in input]
        pd.concat(dfs).to_csv(output.o, index=False, sep='\t')
