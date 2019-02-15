"""Polygenic risk score after imputation"""

# https://www.nature.com/articles/ng.3359#supplementary-information
rule prep_gwas_base_eur:
    input:
        i = DATA + 'raw/ibd_gwas/ng.3359-S4.xlsx'
    output:
        o = DATA + 'interim/ibd_gwas.eur.assoc'
    run:
        def fix_positions(row):
            if row['SNP']=='rs3172494' and row['BP']==48681053:
                return 48731487
            if row['SNP']=='rs9868809' and row['BP']==48731487:
                return 48681053
            if row['SNP']=='rs4768236' and row['BP']==40528432:
                return 40756472
            if row['SNP']=='rs11564258' and row['BP']==40756472:
                return 40792300
            return row['BP']
        df = pd.read_excel(input.i, sheet_name='Heterogeneity of effect', skiprows=7)
        cols = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'EUR_OR', 'EUR_PVAL', 'EUR_SE']
        df.loc[:, 'BP'] = df.apply(fix_positions, axis=1)
        crit = df.apply(lambda row: row['SNP'] != 'rs2226628', axis=1)
        df[crit][cols].rename(columns={'EUR_OR':'OR', 'EUR_PVAL':'PVAL', 'EUR_SE':'SE'}).to_csv(output.o, index=False, sep=' ')

# use eur snps for now
rule prep_gwas_base_tpop:
    input:
        i = DATA + 'raw/ibd_gwas/ng.3359-S4.xlsx'
    output:
        o = DATA + 'interim/ibd_gwas.tpop.assoc'
    run:
        def fix_positions(row):
            if row['SNP']=='rs3172494' and row['BP']==48681053:
                return 48731487
            if row['SNP']=='rs9868809' and row['BP']==48731487:
                return 48681053
            if row['SNP']=='rs4768236' and row['BP']==40528432:
                return 40756472
            if row['SNP']=='rs11564258' and row['BP']==40756472:
                return 40792300
            return row['BP']
        df = pd.read_excel(input.i, sheet_name='Heterogeneity of effect', skiprows=7)
        cols = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'EUR_OR', 'EUR_PVAL', 'EUR_SE']
        df.loc[:, 'BP'] = df.apply(fix_positions, axis=1)
        crit = df.apply(lambda row: row['SNP'] != 'rs2226628', axis=1)
        df[crit][cols].rename(columns={'EUR_OR':'OR', 'EUR_PVAL':'PVAL', 'EUR_SE':'SE'}).to_csv(output.o, index=False, sep=' ')

rule mk_prs_list:
    input:
        i = DATA + 'interim/ibd_gwas.{pop}.assoc',
        ann = DATA + 'interim/tables_tmp/adult.all.{pop}.assoc.paths.csv'
    output:
        o = DATA + 'interim/prs_adult_ls/{go}.{pop}.assoc'
    run:
        if wildcards.go=='GO_ALL':
            shell('cp {input.i} {output}')
        else:
            df = pd.read_csv(input.i, delim_whitespace=True)
            ann_df = pd.read_csv(input.ann, sep='\t')
            crit = ann_df.apply(lambda row: wildcards.go in row['pathways'], axis=1)
            ann_df = ann_df[crit]
            keep_keys = {row['SNP'] + ':' + row['A1_plink'] + ':' + row['A2_plink'] for _, row in ann_df.iterrows()}
            crit = df.apply(lambda row: ':'.join([str(row[x]) for x in ('CHR', 'BP', 'A1', 'A2')]) in keep_keys, axis=1)
            df[crit].to_csv(output.o, index=False, sep='\t')

rule mk_prsice_sample_ls:
    input:
        m = DATA + 'processed/MANIFEST.csv',
        k = DATA + 'interim/mds_cut/{pop}.keep_samples'
    output:
        keep = DATA + 'interim/sample_subsets/{pop}.{group}',
    run:
        keep_samples = {}
        with open(input.k) as f:
            for line in f:
                keep_samples[line.split()[1].strip()] = True

        def filter_samples(row):
            if not row['IID'] in keep_samples:
                return False

            if wildcards.group=='all':
                return True
            if wildcards.group=='ibd_all':
                # must recode fam file
                return row['HC or IBD or ONC'] != 'HC'
            if wildcards.group=='early':
                return row['HC or IBD or ONC'] == 'HC' or 'VEO' in row['Study Group']
            if wildcards.group=='late':
                return not 'VEO' in row['Study Group']

        df = pd.read_csv(input.m)
        crit = df.apply(filter_samples, axis=1)
        df.loc[:, 'iid'] = df.apply(lambda row: '0_' + row['IID'], axis=1)
        df.loc[:, 'junk'] = 0
        df[crit][['junk', 'iid']].to_csv(output.keep, sep=' ', header=None, index=False)

# split into groups and recode fam based on groups
rule mk_prsice_sample_subsets:
    input:
        keep = DATA + 'interim/sample_subsets/{pop}.{group}',
        b = DATA + 'processed/bfiles_imputed/{pop}.fam',
    output:
        b = DATA + 'interim/bfiles_imputed_grouped_tmp/{group}/{pop}.fam'
    singularity:
        PLINK
    log:
        LOG + 'prs/keep_samples.{group}.{pop}'
    shell:
        """plink --bfile $(dirname {input.b})/{wildcards.pop} \
        --keep {input.keep} --make-bed --out $(dirname {output})/{wildcards.pop} &> {log}"""

rule recode_fam_prsice_sample_subsets:
    """Change pheno status for group comparison"""
    input:
        b = DATA + 'interim/bfiles_imputed_grouped_tmp/{group}/{pop}.fam',
        m = DATA + 'processed/MANIFEST.csv',
    output:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/{pop}.fam'
    run:
        shell('cp $(dirname {input.b})/{wildcards.pop}.bim $(dirname {output.b})/{wildcards.pop}.bim')
        shell('cp $(dirname {input.b})/{wildcards.pop}.bed $(dirname {output.b})/{wildcards.pop}.bed')
        if wildcards.group != 'ibd_all':
            shell('cp $(dirname {input.b})/{wildcards.pop}.fam $(dirname {output.b})/{wildcards.pop}.fam')
        else:
            df = pd.read_csv(input.m)
            df.loc[:, 'pheno'] = df.apply(lambda row: 2 if 'VEO' in row['Study Group'] else 1, axis=1)
            df.loc[:, 'iid'] = df.apply(lambda row: '0_' + row['IID'], axis=1)
            phenos = {row['iid']:row['pheno'] for _, row in df.iterrows()}
            with open(input.b) as fam, open(output.b, 'w') as fout:
                for line in fam:
                    sp = line.strip().split()
                    sample = sp[1]
                    pheno = str(phenos[sample])
                    print(' '.join(sp[:-1] + [pheno]), file=fout)

rule base_impute_overlap:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/{pop}.bim',
        a = DATA + 'interim/prs_adult_ls/{go}.{pop}.assoc',
        init = DATA + 'interim/bfiles_{pop}/3groups.bim',
    output:
        o = DATA + 'interim/prsice/snp_overlap/{go}/{group}.{pop}.imputed_all',
        oi = DATA + 'interim/prsice/snp_overlap/{go}/{group}.{pop}.init'
    run:
        shell('cut -f 1,4 {input.b} | sort -u > {output.o}.b')
        shell('cut -f 1,4 {input.init} | sort -u > {output.o}.init')
        shell('cut -f 1,2 -d " " {input.a} | sed "s/ /\t/g" | sort -u > {output.o}.a')
        shell('comm -12 {output.o}.a {output.o}.b > {output.o}')
        shell('comm -12 {output.o} {output.o}.init > {output.oi}')

rule prsice_eur:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/eur.fam',
        a = DATA + 'interim/prs_adult_ls/{go}.eur.assoc'
    output:
        DATA + 'interim/prsice/{group}/{go}/eur.summary',
        DATA + 'interim/prsice/{group}/{go}/eur.best'
    singularity:
        PRSICE
    threads: 16
    log:
        LOG + 'eur.prs.{go}.{group}'
    shell:
        'Rscript /usr/local/bin/PRSice.R --dir {DATA}/interim/prsice/{wildcards.group}/{wildcards.go}/ '
        '--snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --se SE --pvalue PVAL '
        '--prsice /usr/local/bin/PRSice_linux --keep-ambig '
        '--base {input.a} --perm 1000000 --no-clump '
        '--target {DATA}interim/bfiles_imputed_grouped/{wildcards.group}/eur '
        '--thread {threads} --binary-target T '
        '--out {DATA}interim/prsice/{wildcards.group}/{wildcards.go}/eur &> {log} || touch {log}'

rule mk_covfile:
    input:
        pcs = DATA + 'interim/mds_dat/ibd_hapmap.dat'
    output:
        o = DATA + 'interim/prs/tpop.covfile'
    run:
        df = pd.read_csv(input.pcs, sep='\t')
        df.loc[:, 'IID'] = df.apply(lambda row: '0_' + row['IID'], axis=1)
        df[['FID', 'IID', 'C1', 'C2']].to_csv(output.o, index=False, sep='\t')

# use 1st two pcs
rule prsice_tpop:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/tpop.fam',
        a = DATA + 'interim/prs_adult_ls/{go}.tpop.assoc',
        cv = DATA + 'interim/prs/tpop.covfile'
    output:
        DATA + 'interim/prsice/{group}/{go}/tpop.summary',
        DATA + 'interim/prsice/{group}/{go}/tpop.best'
    singularity:
        PRSICE
    threads: 16
    log:
        LOG + 'tpop.prs.{go}.{group}'
    shell:
        'Rscript /usr/local/bin/PRSice.R --dir {DATA}/interim/prsice/{wildcards.group}/{wildcards.go}/ '
        '--snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --se SE --pvalue PVAL '
        '--prsice /usr/local/bin/PRSice_linux --keep-ambig '
        '--base {input.a} --perm 1000000 --no-clump '
        '--target {DATA}interim/bfiles_imputed_grouped/{wildcards.group}/tpop '
        '--thread {threads} --binary-target T '
        '--cov-file {input.cv} --cov-col "@C[1-2]" '
        '--out {DATA}interim/prsice/{wildcards.group}/{wildcards.go}/tpop &> {log} || touch {log}'

rule annotate_prsice_scores:
    input:
        p = DATA + 'interim/prsice/{group}/{go}/{pop}.best',
        f = DATA + 'interim/bfiles_imputed_grouped/{group}/{pop}.fam',
        m = DATA + 'processed/MANIFEST.csv',
        r = DATA + 'interim/snp_groups/{pop}'
    output:
        o = DATA + 'interim/prsice_parsed/{group}/{go}/{pop}.dat'
    run:
        pheno_dict = { 'all': {1:'HC', 2:'IBD'},
                       'early': {1:'HC', 2:'VEO'},
                       'late': {1:'HC', 2:'Late IBD'},
                       'ibd_all': {1:'Late IBD', 2:'VEO'},
        }

        def recode_pheno(row):
            return pheno_dict[wildcards.group][row['pheno']]

        def mk_subject_id(row):
            si = re.findall(r"\d+", row['SSID'].split('CHOP_')[1])[0]
            return int(si)

        m = pd.read_csv(input.m)[['SubjectID', 'IID', 'HC or IBD or ONC', 'Study Group', 'SSID']].rename(columns={'Study Group':'studyGroup'})
        m.loc[:, 'SubjectID'] = m.apply(lambda row: 'Subj ' + str(int(row['SubjectID'])) if str(row['SubjectID'])!='nan' else 'NA', axis=1)
        r = pd.read_csv(input.r, sep='\t')
        m = pd.merge(m, r, on='IID', how='left')
        m.loc[:, 'IID'] = m.apply(lambda row: '0_' + row['IID'], axis=1)
        cols= ['fid', 'IID', 'f', 'm', 'sex', 'pheno']
        int_cols= ['fid', 'f', 'm', 'sex', 'pheno']
        dtype={'fid':int, 'f':int, 'm':int, 'sex':int, 'pheno':int}
        df_dat = pd.read_csv(input.f, sep=' ', header=None, names=cols, dtype=dtype)
        p = pd.read_csv(input.p, delim_whitespace=True)
        df = pd.merge(df_dat, p, on='IID', how='outer')
        quantile_labels = pd.qcut(df['PRS'], 4, labels=["very-low", "low-medium", "high-medium", "high"])
        df.loc[:, 'quantile'] = quantile_labels
        df.loc[:, 'pheno'] = df.apply(recode_pheno, axis=1)
        df.loc[:, 'GO'] = wildcards.go
        pd.merge(df, m, on='IID', how='left').to_csv(output.o, index=False, sep='\t')

rule plot_prs_quantiles:
    input:
        DATA + 'interim/prsice_parsed/{group}/{go}/{pop}.dat'
    output:
        PLOTS + '{group}.{pop}.{go}.prs.quantiles.png'
    run:
        R("""
        require(ggplot2)
        d = read.delim("{input}", header=TRUE, sep='\t')
        p = ggplot(data=d, aes(y=PRS, x=factor(pheno), group=pheno, colour=studyGroup)) + geom_point() +
        theme_bw() + facet_grid(quantile~., scale="free_y") + ylab('') + xlab('')
        ggsave("{output}", p)
        """)

rule plot_prs_dist:
    input:
        DATA + 'interim/prsice_parsed/{group}/{go}/{pop}.dat'
    output:
        PLOTS + '{group}.{pop}.{go}.prs.density.png'
    run:
        R("""
        require(ggplot2)
        d = read.delim("{input}", header=TRUE, sep='\t')
        p = ggplot(data=d) + geom_density(aes(x=PRS, group=pheno, colour=factor(pheno))) +
        theme_bw()
        ggsave("{output}", p)
        """)

rule tmp_yue:
    input:
        expand(DATA + 'interim/prsice_parsed/{group}/{pop}.dat', group='all', pop=('eur', 'tpop'))

rule calc_prs_roc:
    input:
        i = DATA + 'interim/prsice_parsed/{group}/{go}/{pop}.dat'
    output:
        o = DATA + 'interim/prsice_roc/{group}.{pop}.{go}.roc',
        auc = DATA + 'interim/prsice_roc/{group}.{pop}.{go}.auc'
    run:
        df = pd.read_csv(input.i, sep='\t')
        if wildcards.group != 'ibd_all':
            df.loc[:, 'y'] = df.apply(lambda row: 0 if row['pheno']=='HC' else 1, axis=1)
        else:
            df.loc[:, 'y'] = df.apply(lambda row: 1 if row['pheno']=='VEO' else 0, axis=1)

        precision, recall, thresholds = precision_recall_curve(df['y'], df['PRS'], pos_label=1)
        fpr, tpr, thresholds = roc_curve(df['y'], df['PRS'], pos_label=1)
        auc = metrics.auc(fpr, tpr)
        curve = {'fpr':fpr, 'tpr':tpr}
        s = pd.DataFrame(curve, columns=['pre', 'rec', 'fpr', 'tpr'])
        s.to_csv(output.o, index=False, sep='\t')
        with open(output.auc, 'w') as fout:
            print('test\tauc', file=fout)
            print('\t'.join((wildcards.group, str(auc))), file=fout)

rule calc_prs_roc_pval:
    input:
        i = DATA + 'interim/prsice_parsed/{group}/{go}/{pop}.dat'
    output:
        o = DATA + 'interim/prsice_roc/{group}.{pop}.{go}.roc_pval',
    conda:
        ENVS + 'verification-env.yml'
    shell:
        "Rscript {SCRIPTS}/roc_pval.R {input} {output}"

rule plot_prs_roc:
    input:
        DATA + 'interim/prsice_roc/{group}.{pop}.{go}.roc'
    output:
        PLOTS + '{group}.{pop}.{go}.prs.roc.png'
    run:
        R("""require(ggplot2)
             d = read.delim("{input}", sep='\t', header=TRUE)
             p = ggplot(data=d) + geom_line(aes(x=fpr, y=tpr)) +
                 theme_bw(base_size=18) +
                 xlab('False positive rate') + ylab('True positive rate') +
                 theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
            ggsave("{output}", p)
          """)

rule join_rs_dat:
    input:
        prs = DATA + 'interim/prsice/{group}/{go}/{pop}.summary',
        roc = DATA + 'interim/prsice_roc/{group}.{pop}.{go}.auc'
    output:
        o = DATA + 'interim/prsice/{group}/{go}/{pop}.dat'
    run:
        def read_df(afile):
            df = pd.read_csv(afile, sep='\t')
            df['test'] = wildcards.group
            df['pathway'] = wildcards.go
            return df

        df = read_df(input.prs)
        roc = pd.read_csv(input.roc, sep='\t')
        pd.merge(df, roc, on='test', how='left').to_csv(output.o, index=False, sep='\t')

rule combine_prs:
    input:
        expand(DATA + 'interim/prsice/{group}/{go}/{{pop}}.dat', group=G, go=GOS)
    output:
        o = PWD + 'writeup/tables/prs.{pop}.md'
    run:
        cols = ['test', 'pathway', 'Coefficient', 'Standard.Error', 'P', 'Empirical-P', 'Num_SNP', 'auc']
        def read_df(afile):
            df = pd.read_csv(afile, sep='\t')[cols]
            return df
        df = pd.concat([read_df(af) for af in input]).sort_values(by='Empirical-P', ascending=True)
        with open(output.o, 'w') as fout:
            print(tabulate.tabulate(df.values, df.columns, tablefmt="pipe"), file=fout)
