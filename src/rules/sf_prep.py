"""Prep CAG PLINK data"""

def mk_raw_bfiles(wc):
    if wc.group=='44':
        return DATA + 'raw/plink44/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped.fam'
    elif wc.group=='185':
        return DATA + 'raw/plink185/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped.fam'

rule discard_samples:
    """These should not be in study."""
    input:
        f = mk_raw_bfiles,
        d = DATA + 'processed/DISCARD_SAMPLES'
    output:
        DATA + 'interim/bfiles_cag/{group}.fam'
    singularity:
        PLINK
    log:
        LOG + 'prep/discard.{group}'
    shell:
        """plink --bfile $(dirname {input.f})/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped \
        --remove {input.d} --make-bed --out $(dirname {output})/{wildcards.group} &> {log}"""

rule fix_bfiles:
    """Fill in pheno and gender.
       Missing samples.
    """
    input:
        i = DATA + 'interim/bfiles_cag/{group}.fam',
        m = DATA + 'processed/MANIFEST.csv'
    output:
        expand(DATA + 'interim/bfiles/{{group,44|185}}.{suffix}', suffix = ('fam', 'bim', 'bed') ),
        DATA + 'processed/onc_samples.{group}'
    run:
        bim = input.i.replace('.fam', '.bim')
        shell('cp {bim} {DATA}interim/bfiles/{wildcards.group}.bim')

        bed = input.i.replace('.fam', '.bed')
        shell('cp {bed} {DATA}interim/bfiles/{wildcards.group}.bed')

        def calc_sex(row):
            if row['gender']=='Male':
                return '1'
            elif row['gender']=='Female':
                return '2'
            else:
                i = 1/0
                print('no sex', row)
                return 'wtf'

        df = pd.read_csv(input.m)
        df.loc[:, 'sex'] = df.apply(calc_sex, axis=1)
        sex = {row['IID']:row['sex'] for _, row in df.iterrows()}
        pheno = {row['IID']:row['HC or IBD or ONC'] for _, row in df.iterrows()}
        trans = {row['IID']:row['SSID'] for _, row in df.iterrows()}
        with open(input.i) as f, open(DATA + 'interim/bfiles/' + wildcards.group + '.fam', 'w') as fout, open(DATA + 'processed/onc_samples.' + wildcards.group, 'w') as fonc:
            for line in f:
                sp = line.strip().split()
                iid = sp[1]

                if pheno[iid]=='HC':
                    p = '1'
                elif pheno[iid] == 'IBD':
                    p = '2'
                elif pheno[iid] == 'ONC':
                    p = '1'
                    print('onc??', wildcards.group, iid, pheno[iid], trans[iid])
                else:
                    print(iid, pheno[iid])
                    i = 1/0

                # collect samples to remvoe later
                if pheno[iid]=='ONC':
                    print(' '.join(sp[:-2] + [sex[iid], p]), file=fonc)

                print(' '.join(sp[:-2] + [sex[iid], p]), file=fout)

rule rm_dups_and_onc:
    input:
        o1 = DATA + 'processed/onc_samples.44',
        o2 = DATA + 'processed/onc_samples.185',
        d = DATA + 'raw/conrad/dup_samples',
        b = expand(DATA + 'interim/bfiles/{group}.fam', group = (44, 185) ),
    output:
        DATA + 'interim/bfiles_rm_samples/{group}.fam'
    singularity:
        PLINK
    log:
        LOG + 'prep/discard_onc.{group}'
    shell:
        "cat {input.o1} {input.o2} {input.d} > tmp_samples && "
        "plink --bfile  {DATA}interim/bfiles/{wildcards.group} "
        "--remove tmp_samples --make-bed --out $(dirname {output})/{wildcards.group} &> {log} && "
        " rm tmp_samples"

rule combine_bfiles:
    input:
        f = expand(DATA + 'interim/bfiles_rm_samples/{group}.fam', group=(44, 185)),
        m = DATA + 'processed/MANIFEST.csv'
    output:
        DATA + 'interim/bfiles/3groups.fam'
    singularity:
        PLINK
    log:
        LOG + 'prep/combine_44_185'
    shell:
        "plink --bfile {DATA}interim/bfiles_rm_samples/185 --bmerge {DATA}interim/bfiles_rm_samples/44 "
        "--chr 1-22, x --make-bed --out {DATA}interim/bfiles/3groups &> {log}"
