"""Polygenic risk score after imputation"""

# https://www.nature.com/articles/ng.3359#supplementary-information
rule prep_gwas_base:
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
        df[cols].rename(columns={'EUR_OR':'OR', 'EUR_PVAL':'PVAL', 'EUR_SE':'SE'}).to_csv(output.o, index=False, sep=' ')

rule mk_prsice_sample_ls:
    input:
        m = DATA + 'processed/MANIFEST.csv',
        k = DATA + 'interim/mds_cut/3groups.keep_samples'
    output:
        keep = DATA + 'interim/sample_subsets.{group}',
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
        keep = DATA + 'interim/sample_subsets.{group}',
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
        a = DATA + 'interim/ibd_gwas.{pop}.assoc',
        init = DATA + 'interim/bfiles_{pop}/3groups.bim',
    output:
        o = DATA + 'interim/prsice/snp_overlap/{group}.{pop}.imputed_all',
        oi = DATA + 'interim/prsice/snp_overlap/{group}.{pop}.init'
    run:
        shell('cut -f 1,4 {input.b} | sort -u > {output.o}.b')
        shell('cut -f 1,4 {input.init} | sort -u > {output.o}.init')
        shell('cut -f 1,2 -d " " {input.a} | sed "s/ /\t/g" | sort -u > {output.o}.a')
        shell('comm -12 {output.o}.a {output.o}.b > {output.o}')
        shell('comm -12 {output.o} {output.o}.init > {output.oi}')

rule prsice:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/{pop}.fam',
        a = DATA + 'interim/ibd_gwas.{pop}.assoc'
    output:
        DATA + 'interim/prsice/{group}/{pop}.summary',
        DATA + 'interim/prsice/{group}/{pop}.best'
    singularity:
        PRSICE
    threads: 16
    log:
        LOG + '{pop}.prs.{group}'
    shell:
        'Rscript /usr/local/bin/PRSice.R --dir {DATA}/interim/prsice/{wildcards.group}/ '
        '--snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --se SE --pvalue PVAL '
        '--prsice /usr/local/bin/PRSice_linux --keep-ambig '
        '--base {input.a} --perm 1000000 --no-clump '
        '--target {DATA}interim/bfiles_imputed_grouped/{wildcards.group}/{wildcards.pop} '
        '--thread {threads} --binary-target T '
        '--out {DATA}interim/prsice/{wildcards.group}/{wildcards.pop} &> {log}'

rule annotate_prsice_scores:
    input:
        p = DATA + 'interim/prsice/{group}/{pop}.best',
        f = DATA + 'interim/bfiles_imputed_grouped/{group}/{pop}.fam',
        m = DATA + 'processed/MANIFEST.csv'
    output:
        o = DATA + 'interim/prsice_parsed/{group}/{pop}.dat'
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

        m = pd.read_csv(input.m)[['IID', 'HC or IBD or ONC', 'Study Group', 'SSID']].rename(columns={'Study Group':'studyGroup'})
        m.loc[:, 'IID'] = m.apply(lambda row: '0_' + row['IID'], axis=1)
        m.loc[:, 'SubjectID'] = m.apply(mk_subject_id, axis=1)
        cols= ['fid', 'IID', 'f', 'm', 'sex', 'pheno']
        int_cols= ['fid', 'f', 'm', 'sex', 'pheno']
        dtype={'fid':int, 'f':int, 'm':int, 'sex':int, 'pheno':int}
        df_dat = pd.read_csv(input.f, sep=' ', header=None, names=cols, dtype=dtype)
        p = pd.read_csv(input.p, delim_whitespace=True)
        df = pd.merge(df_dat, p, on='IID', how='outer')
        quantile_labels = pd.qcut(df['PRS'], 4, labels=["very-low", "low-medium", "high-medium", "high"])
        df.loc[:, 'quantile'] = quantile_labels
        df.loc[:, 'pheno'] = df.apply(recode_pheno, axis=1)
        pd.merge(df, m, on='IID', how='left').to_csv(output.o, index=False, sep='\t')

rule plot_prs_quantiles:
    input:
        DATA + 'interim/prsice_parsed/{group}/{pop}.dat'
    output:
        PLOTS + '{group}.{pop}.prs.quantiles.png'
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
        DATA + 'interim/prsice_parsed/{group}/{pop}.dat'
    output:
        PLOTS + '{group}.{pop}.prs.density.png'
    run:
        R("""
        require(ggplot2)
        d = read.delim("{input}", header=TRUE, sep='\t')
        p = ggplot(data=d) + geom_density(aes(x=PRS, group=pheno, colour=factor(pheno))) +
        theme_bw()
        ggsave("{output}", p)
        """)

rule calc_prs_roc:
    input:
        i = DATA + 'interim/prsice_parsed/{group}/{pop}.dat'
    output:
        o = DATA + 'interim/prsice_roc/{group}.{pop}.roc',
        auc = DATA + 'interim/prsice_roc/{group}.{pop}.auc'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'y'] = df.apply(lambda row: 1 if row['pheno']==2 else 0, axis=1)
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
        i = DATA + 'interim/prsice_parsed/{group}/{pop}.dat'
    output:
        o = DATA + 'interim/prsice_roc/{group}.{pop}.roc_pval',
    conda:
        ENVS + 'verification-env.yml'
    shell:
        "Rscript {SCRIPTS}/roc_pval.R {input} {output}"

rule plot_prs_roc:
    input:
        DATA + 'interim/prsice_roc/{group}.{pop}.roc'
    output:
        PLOTS + '{group}.{pop}.prs.roc.png'
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
        prs = DATA + 'interim/prsice/{group}/{pop}.summary',
        roc = DATA + 'interim/prsice_roc/{group}.{pop}.auc'
    output:
        o = temp(DATA + 'interim/prsice/{group}/{pop}.dat')
    run:
        def read_df(afile):
            df = pd.read_csv(afile, sep='\t')
            df['test'] = afile.split('/')[-2]
            return df

        df = read_df(input.prs)
        roc = pd.read_csv(input.roc, sep='\t')
        pd.merge(df, roc, on='test', how='left').to_csv(output.o, index=False, sep='\t')

rule combine_prs:
    input:
        expand(DATA + 'interim/prsice/{group}/{{pop}}.dat', group=G)
    output:
        o = PWD + 'writeup/tables/prs.{pop}.md'
    run:
        def read_df(afile):
            df = pd.read_csv(afile, sep='\t')
            return df
        df = pd.concat([read_df(af) for af in input])
        with open(output.o, 'w') as fout:
            print(tabulate.tabulate(df.values, df.columns, tablefmt="pipe"), file=fout)
