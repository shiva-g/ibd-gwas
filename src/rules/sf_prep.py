"""Prep PLINK data"""

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

rule filter_snps:
    """rm snps missing in more than 5% of samples: --geno
       rm non-polymorphic: --min-ac
    """
    input:
        expand(DATA + 'interim/bfiles/3groups.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/bfiles_filter_snps/3groups.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/filter_snps'
    shell:
        "plink --bfile {DATA}interim/bfiles/3groups --min-ac 1 --geno 0.05 "
        "--make-bed --out {DATA}interim/bfiles_filter_snps/3groups &> {log}"

rule filter_samples:
    """rm samples missing in more than 5% of snps: --mind
    """
    input:
        expand(DATA + 'interim/bfiles_filter_snps/3groups.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/bfiles_filter_samples/3groups.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prep/filter_samples'
    shell:
        "plink --bfile {DATA}interim/bfiles/3groups --mind 0.05 --make-bed "
        "--out {DATA}interim/bfiles_filter_samples/3groups &> {log}"

rule missing:
    input:
        expand(DATA + 'interim/bfiles/3groups.{suffix}', suffix=('fam', 'bed', 'bim') )
    output:
        expand(DATA + 'interim/missing/3groups.{miss}', miss=('imiss', 'lmiss') )
    singularity:
        PLINK
    shell:
        "plink --bfile {DATA}interim/bfiles/3groups --missing --out {DATA}interim/missing/3groups"


rule bfiles:
    input: expand(DATA + 'interim/bfiles/{group}.fam', group=(44, 185))
#
#https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
