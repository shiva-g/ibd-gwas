rule combine_hapmap_ibd:
    input:
        expand(DATA + 'interim/bfiles_indep/{group}.fam', group=('3groups', 'hapmap'))
    output:
        DATA + 'interim/bfiles_tmp/hapmap.bim',
        DATA + 'interim/bfiles_tmp/3groupshapmap.bim',
    singularity:
        PLINK
    log:
        LOG + 'prep/combine_ibd_hapmap'
    shell:
        "plink --bfile {DATA}interim/bfiles_indep/3groups --bmerge {DATA}interim/bfiles_indep/hapmap "
        "--out {DATA}interim/bfiles_tmp/ibd_hapmap &> {log}.1 || "
        """grep "Warning: Multiple " {DATA}interim/bfiles_tmp/ibd_hapmap.log | """
        """cut -f 7 -d " " | sed "s/.\{{2\}}$//" | sed "s/^.\{{1\}}//" > {DATA}/interim/bfiles_tmp/multi_pos && """
        "cat {DATA}/interim/bfiles_tmp/multi_pos {DATA}interim/bfiles_tmp/ibd_hapmap.missnp > {DATA}interim/bfiles_tmp/ibd_hapmap.ex && "
        "plink --bfile {DATA}interim/bfiles_indep/3groups --biallelic-only strict --exclude {DATA}interim/bfiles_tmp/ibd_hapmap.ex --make-bed --out {DATA}interim/bfiles_tmp/3groups &> {log}.2 && "
        "plink --bfile {DATA}interim/bfiles_indep/hapmap --biallelic-only strict --exclude {DATA}interim/bfiles_tmp/ibd_hapmap.ex --make-bed --out {DATA}interim/bfiles_tmp/hapmap &> {log}.3 "

rule list_discordant_pos_hapmap_ibd:
    input:
        ibd = DATA + 'interim/bfiles_tmp/3groups.bim',
        hp = DATA + 'interim/bfiles_tmp/hapmap.bim'
    output:
        bout = DATA + 'interim/bfiles_tmp/3groups.discord',
        hout = DATA + 'interim/bfiles_tmp/hapmap.discord',
    run:
        names = ['chrom', 'id', 'blank', 'pos', 'allele1', 'allele2']
        names_hp = ['chrom', 'id_hp', 'blank', 'pos', 'allele1_hp', 'allele2_hp']
        ibd = pd.read_csv(input.ibd, sep='\t', header=None, names=names)
        hp = pd.read_csv(input.hp, sep='\t', header=None, names=names_hp)
        # discord when pos matches, but allele1 and allele2 do not
        m = pd.merge(hp, ibd, on='pos', how='inner')
        m.to_csv('xxx', index=False, sep='\t')
        crit = m.apply(lambda row: row['allele1'] != row['allele1_hp'] or row['allele2'] != row['allele2_hp'], axis=1)
        m[crit][['id']].to_csv(output.bout, index=False, header=None)
        m[crit][['id_hp']].to_csv(output.hout, index=False, header=None)

rule combine_hapmap_ibd_2:
    input:
        b=expand(DATA + 'interim/bfiles_tmp/{group}.fam', group=('3groups', 'hapmap')),
        ibd_rm = DATA + 'interim/bfiles_tmp/3groups.discord',
        hp_rm = DATA + 'interim/bfiles_tmp/hapmap.discord',
    output:
        DATA + 'interim/bfiles/ibd_hapmap.fam',
        DATA + 'interim/bfiles_tmp2/3groups.bim'
    singularity:
        PLINK
    log:
        LOG + 'prep/combine_ibd_hapmap_2'
    shell:
        "plink --bfile {DATA}interim/bfiles_tmp/3groups --exclude {DATA}interim/bfiles_tmp/3groups.discord --make-bed --out {DATA}interim/bfiles_tmp2/3groups &> {log}.4 && "
        "plink --bfile {DATA}interim/bfiles_tmp/hapmap -exclude {DATA}interim/bfiles_tmp/hapmap.discord --make-bed --out {DATA}interim/bfiles_tmp2/hapmap &> {log}.5 && "
        "plink --bfile {DATA}interim/bfiles_tmp2/3groups --bmerge {DATA}interim/bfiles_tmp2/hapmap "
        "--merge-mode 1 --merge-equal-pos --make-bed --out {DATA}interim/bfiles/ibd_hapmap &> {log}"

rule ibd_hapmap_ibd:
    input:
        DATA + 'interim/bfiles_filter_samples/{group}.fam'
    output:
        DATA + 'interim/plink_genome/{group}.genome'
    singularity:
        PLINK
    log:
        LOG + 'mds/genome.{group}'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_samples/{wildcards.group} --genome "
        "--out {DATA}interim/plink_genome/{wildcards.group} &> {log}"

rule ibd_hapmap_mds:
    input:
        f = DATA + 'interim/bfiles_filter_samples/{group}.fam',
        g = DATA + 'interim/plink_genome/{group}.genome'
    output:
        DATA + 'interim/plink_mds/{group}.mds'
    singularity:
        PLINK
    log:
        LOG + 'mds/mds_{group}'
    shell:
        "plink --bfile {DATA}interim/bfiles_filter_samples/{wildcards.group} "
        "--read-genome {DATA}interim/plink_genome/{wildcards.group}.genome "
        "--cluster --mds-plot 2 --out {DATA}interim/plink_mds/{wildcards.group} &> {log}"

rule color_mds:
    input:
        mds = DATA + 'interim/plink_mds/{group}.mds',
        yri = DATA + 'interim/hapmap/YRI.fam',
        ceu = DATA + 'interim/hapmap/CEU.fam',
        asn = DATA + 'interim/hapmap/JPT_CHB.fam',
        ibd = DATA + 'interim/bfiles_filter_samples/3groups.fam',
    output:
        o = DATA + 'interim/mds_dat/{group}.dat'
    run:
        mds = pd.read_csv(input.mds, delim_whitespace=True)
        mds.loc[:, 'FID'] = mds.apply(lambda row: row['FID'].lstrip().strip(), axis=1)
        # fam = pd.read_csv(input.f, header=None, names=['FID', 'IID', 'father', 'mother', 'sex', 'phenotype'], delimiter=' ')
        yri_fam, ceu_fam, asn_fam, ibd_fam = [pd.read_csv(f, header=None, names=['FID', 'IID', 'father', 'mother', 'sex', 'phenotype'], delimiter=' ') for f in (input.yri, input.ceu, input.asn, input.ibd)]
        yri_fam['group'] = 'YRI'
        ibd_fam.loc[:, 'FID'] = mds.apply(lambda row: str(row['FID']), axis=1)
        ibd_fam.loc[:, 'group'] = ibd_fam.apply(lambda row: 'case' if row['phenotype']==2 else 'control', axis=1)
        asn_fam['group'] = 'JPT_CHB'
        ceu_fam['group'] = 'CEU'
        fam = pd.concat([ceu_fam, asn_fam, ibd_fam, yri_fam])
        df = pd.merge(mds, fam, on=['IID', 'FID'], how='left')
        df.to_csv(output.o, index=False, sep='\t')

rule cut_mds:
    """Restrict study samples."""
    input:
        i = DATA + 'interim/mds_dat/ibd_hapmap.dat'
    output:
        o = DATA + 'interim/mds_cut/3groups.keep_samples'
    run:
        df = pd.read_csv(input.i, sep='\t')
        df[(df.C1<-.025) & (df.C2>.025) & ((df.group=='case') | (df.group=='control'))][['FID', 'IID']].to_csv(output.o, index=False, header=None, sep=' ')

rule restrict_ibd_samples:
    input:
        b = DATA + 'interim/bfiles_indep/{group}.fam',
        k = DATA + 'interim/mds_cut/{group}.keep_samples'
    output:
        DATA + 'processed/bfiles/{group}.bim',
    singularity:
        PLINK
    log:
        LOG + 'mds/restrict.{group}'
    shell:
        "plink --bfile {DATA}interim/bfiles_indep/{wildcards.group} --keep {input.k} "
        "--make-bed --out {DATA}processed/bfiles/{wildcards.group} &> {log}"

rule ibd_ibd:
    input:
        DATA + 'processed/bfiles/{group}.fam'
    output:
        DATA + 'interim/plink_genome_ibd/{group}.genome'
    singularity:
        PLINK
    log:
        LOG + 'mds/ibd_genome.{group}'
    shell:
        "plink --bfile {DATA}processed/bfiles/{wildcards.group} --genome "
        "--out {DATA}interim/plink_genome_ibd/{wildcards.group} &> {log}"

rule ibd_mds:
    input:
        f = DATA + 'processed/bfiles/{group}.fam',
        g = DATA + 'interim/plink_genome_ibd/{group}.genome'
    output:
        DATA + 'interim/plink_mds_ibd/{group}.mds'
    singularity:
        PLINK
    log:
        LOG + 'mds/mds_ibd_{group}'
    shell:
        "plink --bfile {DATA}processed/bfiles/{wildcards.group} "
        "--read-genome {DATA}interim/plink_genome_ibd/{wildcards.group}.genome "
        "--cluster --mds-plot 2 --out {DATA}interim/plink_mds_ibd/{wildcards.group} &> {log}"

rule color_mds_ibd:
    input:
        mds = DATA + 'interim/plink_mds_ibd/{group}.mds',
        ibd = DATA + 'interim/bfiles_filter_samples/{group}.fam',
    output:
        o = DATA + 'interim/mds_dat_ibd/{group}.dat'
    run:
        mds = pd.read_csv(input.mds, delim_whitespace=True)
        fam = pd.read_csv(input.ibd, header=None, names=['FID', 'IID', 'father', 'mother', 'sex', 'phenotype'], delimiter=' ')
        fam.loc[:, 'group'] = fam.apply(lambda row: 'case' if row['phenotype']==2 else 'control', axis=1)
        df = pd.merge(mds, fam, on=['IID', 'FID'], how='left')
        df.to_csv(output.o, index=False, sep='\t')

rule mds:
    input: DATA + 'interim/mds_dat_ibd/3groups.dat'
