"""Polygenic risk score after imputation"""

# https://www.nature.com/articles/ng.3359#supplementary-information
rule prep_gwas_base:
    input:
        i = DATA + 'raw/ibd_gwas/ng.3359-S4.xlsx'
    output:
        o = DATA + 'interim/ibd_gwas.assoc'
    run:
        df = pd.read_excel(input.i, sheet_name='Heterogeneity of effect', skiprows=7)
        cols = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'EUR_OR', 'EUR_PVAL', 'EUR_SE']
        df[cols].to_csv(output.o, index=False, sep=' ')

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
        b = DATA + 'processed/bfiles_imputed/eur.fam',
    output:
        b = DATA + 'interim/bfiles_imputed_grouped_tmp/{group}/eur.fam'
    singularity:
        PLINK
    log:
        LOG + 'prs/keep_samples.{group}'
    shell:
        """plink --bfile $(dirname {input.b})/eur \
        --keep {input.keep} --make-bed --out $(dirname {output})/eur &> {log}"""

rule recode_fam_prsice_sample_subsets:
    """Change pheno status for group comparison"""
    input:
        b = DATA + 'interim/bfiles_imputed_grouped_tmp/{group}/eur.fam',
        m = DATA + 'processed/MANIFEST.csv',
    output:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/eur.fam'
    run:
        shell('cp $(dirname {input.b})/eur.bim $(dirname {output.b})/eur.bim')
        shell('cp $(dirname {input.b})/eur.bed $(dirname {output.b})/eur.bed')
        if wildcards.group != 'ibd_all':
            shell('cp $(dirname {input.b})/eur.fam $(dirname {output.b})/eur.fam')
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

rule prsice:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/eur.fam',
        a = DATA + 'interim/ibd_gwas.assoc'
    output:
        DATA + 'interim/prsice/{group}/eur.summary',
        DATA + 'interim/prsice/{group}/eur.best'
    singularity:
        PRSICE
    threads: 16
    log:
        LOG + 'eur.prs.{group}'
    shell:
        'Rscript /usr/local/bin/PRSice.R --dir {DATA}/interim/prsice/{wildcards.group}/ '
        '--snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat EUR_OR --se EUR_SE --pvalue EUR_PVAL '
        '--prsice /usr/local/bin/PRSice_linux --keep-ambig '
        '--base {input.a} --perm 1000000 --no-clump '
        '--target {DATA}interim/bfiles_imputed_grouped/{wildcards.group}/eur '
        '--thread {threads} --binary-target T '
        '--out {DATA}interim/prsice/{wildcards.group}/eur &> {log}'

rule annotate_prsice_scores:
    input:
        p = DATA + 'interim/prsice/{group}/eur.best',
        f = DATA + 'interim/bfiles_imputed_grouped/{group}/eur.fam'
    output:
        o = DATA + 'interim/prsice_parsed/{group}/eur.dat'
    run:
        cols= ['fid', 'IID', 'f', 'm', 'sex', 'pheno']
        int_cols= ['fid', 'f', 'm', 'sex', 'pheno']
        dtype={'fid':int, 'f':int, 'm':int, 'sex':int, 'pheno':int}
        df_dat = pd.read_csv(input.f, sep=' ', header=None, names=cols, dtype=dtype)
        p = pd.read_csv(input.p, delim_whitespace=True)
        df = pd.merge(df_dat, p, on='IID', how='outer')
        df.to_csv(output.o, index=False, sep='\t')

rule plot_prs_dist:
    input:
        DATA + 'interim/prsice_parsed/{group}/eur.dat'
    output:
        PLOTS + '{group}.eur.prs.density.png'
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
        i = DATA + 'interim/prsice_parsed/{group}/eur.dat'
    output:
        o = DATA + 'interi/prsice_roc/{group}.eur.roc'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df.loc[:, 'y'] = df.apply(lambda row: 1 if row['pheno']==2 else 0, axis=1)
        precision, recall, thresholds = precision_recall_curve(df['y'], df['PRS'], pos_label=1)
        fpr, tpr, thresholds = roc_curve(df['y'], df['PRS'], pos_label=1)
        curve = {'fpr':fpr, 'tpr':tpr}
        s = pd.DataFrame(curve, columns=['pre', 'rec', 'fpr', 'tpr'])
        s.to_csv(output.o, index=False, sep='\t')

rule plot_prs_roc:
    input:
        DATA + 'interi/prsice_roc/{group}.eur.roc'
    output:
        PLOTS + '{group}.eur.prs.roc.png'
    run:
        R("""require(ggplot2)
             d = read.delim("{input}", sep='\t', header=TRUE)
             p = ggplot(data=d) + geom_line(aes(x=fpr, y=tpr)) +
                 theme_bw(base_size=18) +
                 xlab('False positive rate') + ylab('True positive rate') +
                 theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
            ggsave("{output}", p)
          """)

rule combine_prs:
    input:
        expand(DATA + 'interim/prsice/{group}/eur.summary', group=G)
    output:
        o = PWD + 'writeup/tables/prs.md'
    run:
        def read_df(afile):
            df = pd.read_csv(afile, sep='\t')
            df['test'] = afile.split('/')[-2]
            return df
        df = pd.concat([read_df(af) for af in input])
        with open(output.o, 'w') as fout:
            print(tabulate.tabulate(df.values, df.columns, tablefmt="pipe"), file=fout)
