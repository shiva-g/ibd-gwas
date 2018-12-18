"""Make imputed vcf files into bfiles w/ correct snp names and fam data.
   Only use good imputations (r2>0.3) and maf>1%
"""

rule unzip_imputed:
    input:
        DATA + 'interim/imputed_vcf/chr_{c}.zip'
    output:
        DATA + 'interim/imputed_vcf/chr{c}.dose.vcf.gz'
    shell:
        'unzip -o -P "Ok4iEYdK9Ue\Pu" -d {DATA}/interim/imputed_vcf/ {input}'

rule limit_imputed_r2:
    input:
        DATA + 'interim/imputed_vcf/chr{c}.dose.vcf.gz'
    output:
        DATA + 'interim/imputed_r2_limit_vcf/chr{c}.vcf.gz'
    conda:
        ENVS + 'bcftools-env.yml'
    shell:
        'bcftools filter --include "INFO/R2>0.3" {input} | bgzip -c > {output}'

rule vcf_to_bfiles:
    input:
        DATA + 'interim/imputed_r2_limit_vcf/chr{chr}.vcf.gz'
    output:
        expand(DATA + 'interim/bfiles_imputed/chr{{chr}}.{suffix}', suffix=('fam', 'bed', 'bim') )
    singularity:
        PLINK
    log:
        LOG + 'prs/{chr}.vcf'
    shell:
        "plink --vcf {input} --make-bed --const-fid 0 "
        "--out {DATA}interim/bfiles_imputed/chr{wildcards.chr} &> {log}"

# rule imputed_maf_cutoff:
#     input:
#         expand(DATA + 'interim/bfiles_imputed_tmp/chr{{chr}}.{suffix}', suffix=('fam', 'bed', 'bim') ),
#         b = DATA + 'interim/bfiles_imputed_tmp/chr{chr}.bim'
#     output:
#         expand(DATA + 'interim/bfiles_imputed/chr{{chr}}.{suffix}', suffix=('fam', 'bed', 'bim') )
#     singularity:
#         PLINK
#     log:
#         LOG + 'prs/{chr}.vcf'
#     shell:
#         "plink -bfile $(dirname {input.b})/chr{wildcards.chr} "
#         "--make-bed "
#         "--out {DATA}interim/bfiles_imputed/chr{wildcards.chr} &> {log}"

rule merge_imputed_bfiles_ls:
    input:
        expand(DATA + 'interim/bfiles_imputed_ex/chr{chr}.fam', chr=range(1, 23))
    output:
        o = DATA + 'interim/bfiles_imputed.eur.ls'
    run:
        with open(output.o, 'w') as fout:
            for afile in list(input):
                print(afile.split('.')[0], file=fout)

rule merge_imputed_bfiles_ls_1:
    input:
        expand(DATA + 'interim/bfiles_imputed/chr{chr}.fam', chr=range(1, 23))
    output:
        o = DATA + 'interim/bfiles_imputed.eur.ls1'
    run:
        with open(output.o, 'w') as fout:
            for afile in list(input):
                print(afile.split('.')[0], file=fout)

rule merge_imputed_bfiles_tmp:
    input:
        DATA + 'interim/bfiles_imputed.eur.ls1'
    output:
        DATA + 'interim/bfiles_imputed_tmp/eur-merge.missnp'
    singularity:
        PLINK
    log:
        LOG + 'prs/merge_tmp'
    shell:
        "plink --merge-list {input} --make-bed "
        "--out {DATA}interim/bfiles_imputed_tmp/eur &> {log} || touch {output}"

rule rm_imputed_tri_allelic:
    input:
        f = DATA + 'interim/bfiles_imputed/chr{chr}.fam',
        ex = DATA + 'interim/bfiles_imputed_tmp/eur-merge.missnp'
    output:
        DATA + 'interim/bfiles_imputed_ex/chr{chr}.fam'
    singularity:
        PLINK
    log:
        LOG + 'prs/{chr}.rm_tri_tmp'
    shell:
        "plink --bfile {DATA}/interim/bfiles_imputed/chr{wildcards.chr} "
        "--exclude {input.ex} --make-bed "
        "--out {DATA}interim/bfiles_imputed_ex/chr{wildcards.chr} &> {log} "

rule merge_imputed_bfiles:
    input:
        DATA + 'interim/bfiles_imputed.eur.ls'
    output:
        DATA + 'interim/bfiles_imputed_combined/eur.fam'
    singularity:
        PLINK
    log:
        LOG + 'prs/merge'
    shell:
        "plink --merge-list {input} --make-bed "
        "--out {DATA}interim/bfiles_imputed_combined/eur &> {log} "

rule mk_rs_names:
    input:
        i = DATA + 'interim/ibd_gwas.assoc'
    output:
        o = DATA + 'interim/rs_names'
    run:
        df = pd.read_csv(input.i, sep=' ')
        df.loc[:, 'id'] = df.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        df[['id', 'SNP']].to_csv(output.o, index=False, sep=' ', header=None)

rule update_fam:
    """Mk fam for imputed files b/c imputation removed data.
       Some data will be missing for samples removed from study.
       These await a new imputation run.
    """
    input:
        imputed_fam = DATA + 'interim/bfiles_imputed_combined/eur.fam',
        data_fam = DATA + 'interim/bfiles_eur/3groups.fam'
    output:
        o = DATA + 'interim/tmp.fam'
    run:
        cols= ['fid', 'iid', 'f', 'm', 'sex', 'pheno']
        int_cols= ['fid', 'f', 'm', 'sex', 'pheno']
        dtype={'fid':int, 'f':int, 'm':int, 'sex':int, 'pheno':int}
        df_imputed = pd.read_csv(input.imputed_fam, sep=' ', header=None, names=cols, dtype=dtype)[['iid']]
        df_dat = pd.read_csv(input.data_fam, sep=' ', header=None, names=cols, dtype=dtype)
        df_dat.loc[:, 'iid'] = df_dat.apply(lambda row: '0_' + row['iid'], axis=1)
        bad_samples = set(df_imputed['iid']) - set(df_dat['iid'])
        bad_dat = []
        for bad_sample in bad_samples:
            bad_dat.append( [0, bad_sample, 0, 0, -9, -9])
        bad_df = pd.DataFrame(bad_dat, columns=cols)
        df = pd.merge(df_imputed, pd.concat([df_dat, bad_df]), how='left', on='iid')
        df[cols].to_csv(output.o, sep=' ', index=False, header=False)

rule rename_rs_imputed_bfiles:
    """Apply 1% maf cutoff."""
    input:
        fam = DATA + 'interim/tmp.fam',
        f = DATA + 'interim/bfiles_imputed_combined/eur.fam',
        rs_names = DATA + 'interim/rs_names'
    output:
        DATA + 'processed/bfiles_imputed/eur.fam'
    singularity:
        PLINK
    log:
        LOG + 'prs/rename_rs'
    shell:
        "plink --bfile {DATA}/interim/bfiles_imputed_combined/eur "
        "--update-name {input.rs_names} --make-bed --maf 0.01 "
        "--out {DATA}processed/bfiles_imputed/eur &> {log} && "
        """cp {input.fam} {output}"""
