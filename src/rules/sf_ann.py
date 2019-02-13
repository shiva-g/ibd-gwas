"""Annotate imputed vcf file"""

rule vcf_imputed:
    input:
        DATA + 'processed/bfiles_imputed/{group}.fam',
    output:
        DATA + 'interim/bfiles_clean_vcf/{group}.vcf'
    singularity:
        PLINK
    log:
        LOG + 'ann/vcf.{group}'
    shell:
        "plink --bfile {DATA}processed/bfiles_imputed/{wildcards.group} --recode vcf "
        "--out {DATA}interim/bfiles_clean_vcf/{wildcards.group} &> {log}"

rule snpeff:
    input:
        DATA + "interim/bfiles_clean_vcf/{group}.vcf"
    output:
        vcf=DATA + "interim/variants/snpeff/{group}.vcf"
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="GRCh37.75", # reference name (from `snpeff databases`)
        extra="-Xmx32g -Xms16g"  # optional parameters (e.g., max memory 4g)
    singularity:
        "docker://quay.research.chop.edu/evansj/snpeff-docker:grch37"
    conda:
        ENVS + "project.yml"
    threads: 30
    #
    #shell:
    #    "snpEff -noStats -t {threads} -Xmx32g -Xms16g {params.reference} {input} > {output} 2> {log}"
    #
    wrapper:
        "file:///mnt/isilon/cbmi/variome/perry/projects/ext/snakemake-wrappers/bio/snpeff/wrapper.py"

#java -jar SnpSift.jar geneSets -v db/msigDb/msigdb.v3.1.symbols.gmt test.ann.vcf > test.eff.geneSets.vcf
rule parse_vcf_ann:
    input:
        i = DATA + "interim/variants/snpeff/{pop}.vcf"
    output:
        o = DATA + "interim/variants/snpeff_parsed/{pop}"
    run:
        with open(input.i) as f, open(output.o, 'w') as fout:
            print('SNP\tA1_vcf\tA2_vcf\teff', file=fout)
            for line in f:
                if line[0] != '#':
                    sp = line.strip().split('\t')
                    chro, bp = sp[:2]
                    key = chro + ':' + bp
                    ref, alt = sp[3:5]
                    eff = sp[7]
                    print('\t'.join( (key, ref, alt, eff)), file=fout)

rule ann_plink_adult_assoc:
    input:
        a = DATA + 'interim/plink_assoc_fmt/{group}/{pop}.assoc',
        st = DATA + 'interim/{pop}_snptest_final/snptest.out',
        g = DATA + 'interim/ibd_gwas.{pop}.assoc',
        s = DATA + "interim/variants/snpeff_parsed/{pop}"
    output:
        o = DATA + 'interim/adult_assoc_ann/{group}/{pop}.assoc'
    run:
        dat = pd.read_csv(input.s, sep='\t')
        st_cols = {'frequentist_add_pvalue':'P_snptest', 'alleleA':'A1_snptest', 'alleleB':'A2_snptest', 'all_OR':'OR_snptest'}
        snptest = pd.read_csv(input.st, sep='\t')[['frequentist_add_pvalue', 'chromosome', 'position', 'alleleA', 'alleleB', 'cases_maf', 'controls_maf', 'all_OR']].rename(columns=st_cols).fillna(1)
        snptest.loc[:, 'SNP'] = snptest.apply(lambda row: str(row['chromosome']) + ':' + str(row['position']), axis=1)
        assoc = pd.read_csv(input.a, sep='\t').rename(columns={'A1':'A1_plink', 'A2':'A2_plink', 'P':'P_plink', 'OR':'OR_plink'})
        assoc.loc[:, 'SNP'] = assoc.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        full_assoc = pd.merge(assoc, snptest, on='SNP', how='left')
        adult_gwas = pd.read_csv(input.g, delim_whitespace=True).rename(columns={'SNP':'rs', 'A1':'A1_adult', 'A2':'A2_adult', 'OR':'OR_adult', 'PVAL':'PVAL_adult', 'SE':'SE_adult'})
        adult_gwas.loc[:, 'SNP'] = adult_gwas.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        m1 = pd.merge(adult_gwas, full_assoc, on='SNP', how='left').fillna('')
        pd.merge(m1, dat, how='left', on='SNP').to_csv(output.o, index=False, sep='\t')

rule ann_plink_ped_assoc:
    input:
        a = DATA + 'interim/plink_assoc_fmt/{group}/{pop}.assoc',
        st = DATA + 'interim/{pop}_snptest_final/snptest.out',
        g = DATA + 'interim/ibd_gwas.{pop}.assoc',
        s = DATA + "interim/variants/snpeff_parsed/{pop}",
        adult_table = PWD + 'writeup/tables/adult.all.{pop}.assoc.csv'
    output:
        o = DATA + 'interim/ped_assoc_ann/{group}/{pop}.assoc'
    run:
        adult_genes = {x:True for x in pd.read_csv(input.adult_table)['gene'].values if x}
        dat = pd.read_csv(input.s, sep='\t')
        st_cols = {'frequentist_add_pvalue':'P_snptest', 'alleleA':'A1_snptest', 'alleleB':'A2_snptest', 'all_OR':'OR_snptest'}
        snptest = pd.read_csv(input.st, sep='\t')[['frequentist_add_pvalue', 'chromosome', 'position', 'alleleA', 'alleleB', 'cases_maf', 'controls_maf', 'all_OR']].rename(columns=st_cols).fillna(1)
        snptest.loc[:, 'SNP'] = snptest.apply(lambda row: str(row['chromosome']) + ':' + str(row['position']), axis=1)
        assoc = pd.read_csv(input.a, sep='\t').rename(columns={'A1':'A1_plink', 'A2':'A2_plink', 'P':'P_plink', 'OR':'OR_plink'})
        assoc.loc[:, 'SNP'] = assoc.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        full_assoc = pd.merge(assoc, snptest, on='SNP', how='left')
        adult_gwas = pd.read_csv(input.g, delim_whitespace=True).rename(columns={'SNP':'rs', 'A1':'A1_adult', 'A2':'A2_adult', 'OR':'OR_adult', 'PVAL':'PVAL_adult', 'SE':'SE_adult'})
        adult_gwas.loc[:, 'SNP'] = adult_gwas.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        m1 = pd.merge(full_assoc, adult_gwas, on='SNP', how='left').fillna('')
        gg = pd.merge(m1, dat, how='left', on='SNP')
        gg.loc[:, 'gene_tmp'] = gg.apply(lambda row: str(row['eff']).split('ANN=')[1].split('|')[3] if 'ANN=' in str(row['eff']) else 'none', axis=1)
        gg.loc[:, 'in_adult_gene'] = gg.apply(lambda row: row['gene_tmp'] in adult_genes, axis=1)
        crit = gg.apply(lambda row: row['P_snptest']<1e-5 or row['P_plink']<1e-5 or row['in_adult_gene'], axis=1)
        gg[crit].to_csv(output.o, index=False, sep='\t')

rule gene_assoc:
    input:
        a = DATA + 'interim/{age}_assoc_ann/{group}/{pop}.assoc'
    output:
        o = PWD + 'writeup/tables/{age}.{group}.{pop}.assoc.csv'
    run:
        df = pd.read_csv(input.a, sep='\t')
        df.loc[:, 'gene'] = df.apply(lambda row: str(row['eff']).split('ANN=')[1].split('|')[3] if 'ANN=' in str(row['eff']) else 'none', axis=1)
        cols = ['SNP', 'gene', 'A1_plink', 'A2_plink', 'P_plink', 'OR_plink',
                'A1_snptest', 'A2_snptest', 'P_snptest', 'OR_snptest',
                'cases_maf', 'controls_maf',
                'A1_adult', 'A2_adult', 'PVAL_adult', 'OR_adult',
                'A1_vcf', 'A2_vcf', 'eff',]# 'in_adult_gene']
        if 'ped'==wildcards.age:
            cols.append('in_adult_gene')
        df[cols].sort_values(by='P_snptest', ascending=True).to_csv(output.o, index=False, sep=',')
