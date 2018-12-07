"""Prep CAG PLINK data"""

def mk_raw_bfiles(wc):
    if wc.group=='44':
        return DATA + 'raw/plink44/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped.fam'
    elif wc.group=='185':
        return DATA + 'raw/plink185/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped.fam'

rule fix_bfiles:
    """Fill in pheno and gender.
       Missing samples.
    """
    input:  i = mk_raw_bfiles,
            m = DATA + 'processed/MANIFEST.csv'
    output: expand(DATA + 'interim/bfiles/{{group,44|185}}.{suffix}', suffix = ('fam', 'bim', 'bed') )
    run:
        bim = input.i.replace('.fam', '_Forward.bim')
        shell('cp {bim} {DATA}interim/bfiles/{wildcards.group}.bim')

        bed = input.i.replace('.fam', '.bed')
        shell('cp {bed} {DATA}interim/bfiles/{wildcards.group}.bed')

        # if wildcards.group=='44':
        #     pheno = '2' # case
        # elif wildcards.group=='185':
        #     pheno = '1' # control
        # else:
        #     pheno = 1/0
        #sex = '1' # male
        df = pd.read_csv(input.m)
        def calc_sex(row):
            if row['gender']=='Male':
                return '1'
            elif row['gender']=='Female':
                return '2'
            else:
                print('no sex', row)
                return 'wtf'

        df.loc[:, 'sex'] = df.apply(calc_sex, axis=1)
        sex = {row['IID']:row['sex'] for _, row in df.iterrows()}
        pheno = {row['IID']:row['HC or IBD or ONC'] for _, row in df.iterrows()}
        trans = {row['IID']:row['SSID'] for _, row in df.iterrows()}
        with open(input.i) as f, open(DATA + 'interim/bfiles/' + wildcards.group + '.fam', 'w') as fout:
            for line in f:
                sp = line.strip().split()
                iid = sp[1]

                if sex[iid] != 'wtf':
                    if pheno[iid]=='HC':
                        p = '1'
                    elif pheno[iid] == 'IBD':
                        p = '2'
                    elif pheno[iid] == 'ONC':
                        p = '1'
                        print('onc??', wildcards.group, iid, pheno[iid], trans[iid])
                    else:
                        print(iid, pheno[iid])
                        #i = 1/0

                    print(' '.join(sp[:-2] + [sex[iid], p]), file=fout)
                else:
                    print(' '.join(sp[:-2] + ['1', '1']), file=fout)
                    print('miss', iid, trans[iid])

rule combine_bfiles:
    input:
        expand(DATA + 'interim/bfiles/{group}.fam', group=(44, 185))
    output:
        DATA + 'interim/bfiles/3groups.fam'
    singularity:
        PLINK
    log:
        LOG + 'prep/combine_44_185'
    shell:
        "plink --bfile {DATA}interim/bfiles/185 --bmerge {DATA}interim/bfiles/44 "
        "--remove {CONFIG}dup_samples --chr 1-22, x --make-bed --out {DATA}interim/bfiles/3groups &> {log}"
