"""Prep CAG GSA VEO data.
   https://mail.google.com/mail/u/0/#inbox/FMfcgxwBVDKqJbCSJrRCzNkLHvWCjNGc
   data/raw/veo-ibd/
"""

rule fix_veo_gsa_bfiles:
    """Fill in pheno and gender.
    """
    input:
        i = DATA + 'raw/veo-ibd/IJUK_VEOIBDx266_11-10-17PLINKFILES/IJUK_VEOIBDx266_11-10-17.fam',
        m = DATA + 'processed/MANIFEST.csv',
    output:
        expand(DATA + 'interim/gsa_bfiles_fixed/veo-gsa.{suffix}', suffix = ('fam', 'bim', 'bed') )
    run:
        bim = input.i.replace('.fam', '_Forward.bim')
        shell('cp {bim} {DATA}interim/gsa_bfiles_fixed/veo-gsa.bim')

        bed = input.i.replace('.fam', '.bed')
        shell('cp {bed} {DATA}interim/gsa_bfiles_fixed/veo-gsa.bed')

        pheno = '2' # case
        df = pd.read_csv(input.m)
        sexes = {row['IID']:row['gender'] for _, row in df.iterrows()}

        def calc_sex(row):
            if row=='Male':
                return '1'
            elif row=='Female':
                return '2'
            else:
                i = 1/0
                print('no sex', row)

        with open(input.i) as f, open(DATA + 'interim/gsa_bfiles_fixed/veo-gsa.fam', 'w') as fout:
            for line in f:
                sp = line.strip().split('\t')
                s = sexes[sp[1]]
                sex = calc_sex(s)
                print(' '.join(sp[:-2] + [sex, pheno]), file=fout)

rule fix_hc_gsa_bfiles:
    """Fill in pheno.
       No sex info available.
    """
    input:  i = DATA + 'raw/veo-ibd/IJUK_GSA_719_Controls_12-1-17_PlinkFiles/IJUK_GSA_719_Controls_12-1-17.fam'
    output: expand(DATA + 'interim/gsa_bfiles_fixed/hc-gsa.{suffix}', suffix = ('fam', 'bim', 'bed') )
    run:
        bim = input.i.replace('.fam', '_Forward.bim')
        shell('cp {bim} {DATA}interim/gsa_bfiles_fixed/hc-gsa.bim')

        bed = input.i.replace('.fam', '.bed')
        shell('cp {bed} {DATA}interim/gsa_bfiles_fixed/hc-gsa.bed')

        pheno = '1' # control
        sex = '-9'
        with open(input.i) as f, open(DATA + 'interim/gsa_bfiles_fixed/hc-gsa.fam', 'w') as fout:
            for line in f:
                sp = line.strip().split()
                print(' '.join(sp[:-2] + [sex, pheno]), file=fout)

rule discard_hc_gsa_samples:
    """These control samples should not be in study
       b/c they have IBD related phenotypes."""
    input:
        f = DATA + 'interim/gsa_bfiles_fixed/hc-gsa.fam',
        d = DATA + 'processed/DISCARD_SAMPLES.gsa'
    output:
        DATA + 'interim/gsa_bfiles_fixed2/hc-gsa.fam'
    singularity:
        PLINK
    log:
        LOG + 'prep/discard.gsa.hc'
    shell:
        """plink --bfile $(dirname {input.f})/hc-gsa \
        --remove {input.d} --make-bed \
        --out $(dirname {output})/hc-gsa &> {log}"""

rule fake_discard_veo_gsa_samples:
    input:
        f = DATA + 'interim/gsa_bfiles_fixed/veo-gsa.fam',
    output:
        o = DATA + 'interim/gsa_bfiles_fixed2/veo-gsa.fam'
    run:
        i, o = input.f.replace('.fam', '.bim'), output.o.replace('.fam', '.bim')
        shell('cp {i} {o}')
        i, o = input.f.replace('.fam', '.bed'), output.o.replace('.fam', '.bed')
        shell('cp {i} {o}')
        shell("cp {input} {output}")

rule list_discordant_pos_gsa:
    input:
        ibd = DATA + 'interim/gsa_bfiles_fixed2/veo-gsa.bim',
        hc = DATA + 'interim/gsa_bfiles_fixed2/hc-gsa.bim'
    output:
        bout = DATA + 'interim/gsa_discord/veo.discord',
        hout = DATA + 'interim/gsa_discord/hc.discord',
    run:
        names = ['chrom', 'id', 'blank', 'pos', 'allele1', 'allele2']
        names_hp = ['chrom', 'id_hp', 'blank', 'pos', 'allele1_hp', 'allele2_hp']
        ibd = pd.read_csv(input.ibd, sep='\t', header=None, names=names)
        hp = pd.read_csv(input.hc, sep='\t', header=None, names=names_hp)
        # discord when pos matches, but allele1 and allele2 do not
        m = pd.merge(hp, ibd, on='pos', how='inner')
        #m.to_csv('xxx', index=False, sep='\t')
        crit = m.apply(lambda row: row['allele1'] != row['allele1_hp'] or row['allele2'] != row['allele2_hp'], axis=1)
        m[crit][['id']].to_csv(output.bout, index=False, header=None)
        m[crit][['id_hp']].to_csv(output.hout, index=False, header=None)

rule combine_gsa_cag_bfiles:
    input:
        expand(DATA + 'interim/gsa_bfiles_fixed2/{group}.fam', group=('hc-gsa', 'veo-gsa')),
        ex_hc = DATA + 'interim/gsa_discord/hc.discord',
        ex_veo = DATA + 'interim/gsa_discord/veo.discord',
    output:
        f=DATA + 'interim/bfiles/gsa.fam',
        a=expand(DATA + 'interim/bfiles/gsa.{s}', s=('bim', 'bed',))
    singularity:
        PLINK
    log:
        LOG + 'prep/combine_gsa'
    shell:
        "plink --bfile {DATA}interim/gsa_bfiles_fixed2/veo-gsa "
        "--exclude {input.ex_veo} "
        "--make-bed --out {DATA}interim/bfiles/veo-gsa &> {log}.v && "
        "plink --bfile {DATA}interim/gsa_bfiles_fixed2/hc-gsa -exclude {input.ex_hc} "
        "--make-bed --out {DATA}interim/bfiles/hc-gsa &> {log}.h && "
        "plink --bfile {DATA}interim/bfiles/veo-gsa "
        "--bmerge {DATA}interim/bfiles/hc-gsa --merge-mode 1 --merge-equal-pos "
        "--allow-no-sex --chr 1-22, x --make-bed --out $(dirname {output.f})/gsa &> {log}"
