"""Prep CAG GSA VEO data.
   https://mail.google.com/mail/u/0/#inbox/FMfcgxwBVDKqJbCSJrRCzNkLHvWCjNGc
   data/raw/veo-ibd/
"""
def mk_raw_bfiles(wc):
    if wc.group=='44':
        return DATA + 'raw/plink44/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped.fam'
    elif wc.group=='185':
        return DATA + 'raw/plink185/paths.txt.GSA-mGluR_enrichd_20011739X343186_B2.ped.fam'

 rule fix_bfiles:
    """Fill in pheno.
       Make everyone male
    """
    input:  i = mk_raw_bfiles
    output: expand(DATA + 'interim/bfiles/{{group,44|185}}.{suffix}', suffix = ('fam', 'bim', 'bed') )
    run:
        bim = input.i.replace('.fam', '_Forward.bim')
        shell('cp {bim} {DATA}interim/bfiles/{wildcards.group}.bim')

         bed = input.i.replace('.fam', '.bed')
        shell('cp {bed} {DATA}interim/bfiles/{wildcards.group}.bed')

         if wildcards.group=='44':
            pheno = '2' # case
        elif wildcards.group=='185':
            pheno = '1' # control
        else:
            pheno = 1/0
        sex = '1' # male
        with open(input.i) as f, open(DATA + 'interim/bfiles/' + wildcards.group + '.fam', 'w') as fout:
            for line in f:
                sp = line.strip().split()
                print(' '.join(sp[:-2] + [sex, pheno]), file=fout)

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
        "--autosome-xy --make-bed --out {DATA}interim/bfiles/3groups &> {log}"

rule fix_gsa_cag_bfiles:
    """Fill in pheno and gender.
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

