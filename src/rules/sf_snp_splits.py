"""Split samples by SNP.
   For Yue and Maire.
"""

rule mk_impt_snps:
    input:
        s = DATA + 'raw/important_snps',
        b = DATA + 'processed/bfiles_imputed/{pop}.bim',
    output:
        o = DATA + 'interim/{pop}.snps_to_extract'
    run:
        keys = {}
        with open(input.s) as f:
            f.readline()
            for line in f:
                keys[line.split(',')[0]] = True
        with open(input.b) as f, open(output.o, 'w') as fout:
            for line in f:
                chrom, rs, junk, pos, a1, a2 = line.strip().split('\t')
                key = chrom + ':' + pos
                if key in keys:
                    print(rs, file=fout)

rule extract_snps:
    input:
        f = DATA + 'processed/bfiles_imputed/{pop}.fam',
        s = DATA + 'interim/{pop}.snps_to_extract'
    output:
        p = DATA + 'interim/extracted_snps/{pop}.ped',
        m = DATA + 'interim/extracted_snps/{pop}.map'
    singularity:
        PLINK
    log:
        LOG + 'extracted_snps/{pop}'
    shell:
        'plink --bfile $(dirname {input.f})/{wildcards.pop} '
        '--extract {input.s} --recode '
        '--out $(dirname {output.p})/{wildcards.pop} &> {log}'

rule recode_snp_splits:
    input:
        p = DATA + 'interim/extracted_snps/{pop}.ped',
        m = DATA + 'interim/extracted_snps/{pop}.map',
        s = DATA + 'raw/important_snps',
        man = DATA + 'processed/MANIFEST.csv'
    output:
        o = DATA + 'interim/snp_groups/{pop}'
    run:
        snp_dat = pd.read_csv(input.s)
        risk_alleles = { row['gene']:row['risk_allele'].split('/') for _, row in snp_dat.iterrows() }
        map_df = pd.read_csv(input.m, header=None, sep='\t', names=['chrom', 'rs', 'j', 'pos'])
        map_df.loc[:, 'chrpos'] = map_df.apply(lambda row: str(row['chrom']) + ':' + str(row['pos']), axis=1)
        df = pd.merge(map_df, snp_dat, on='chrpos', how='left')
        cols = ['j', 'id', 'j1', 'j2', 'j3', 'j4']
        genes = df['gene'].values
        for gene in genes:
            cols += [gene + '_a1', gene + '_a2']
        df = pd.read_csv(input.p, header=None, delim_whitespace=True, names=cols)
        df.loc[:, 'IID'] = df.apply(lambda row: row['id'][2:], axis=1)

        def score_risk(row, gene):
            alleles = risk_alleles[gene]
            a1, a2 = row[gene + '_a1'], row[gene + '_a2']
            if a1 in alleles and a2 in alleles:
                return 'hom-risk'
            if a1 in alleles or a2 in alleles:
                return 'het-risk'
            return 'no-risk'

        for gene in genes:
            df.loc[:, gene + '_risk'] = df.apply(lambda row: score_risk(row, gene), axis=1)
        cols = ['IID'] + [gene + '_risk' for gene in genes]
        df[cols].to_csv(output.o, index=False, sep='\t')
