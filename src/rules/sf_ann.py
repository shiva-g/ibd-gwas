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

rule ann_plink_assoc:
    input:
        a = DATA + 'interim/plink_assoc_fmt/{group}/{pop}.assoc',
        g = DATA + 'interim/ibd_gwas.{pop}.assoc',
        s = DATA + "interim/variants/snpeff_parsed/{pop}"
    output:
        o = DATA + 'interim/plink_assoc_fmt_ann/{group}/{pop}.assoc'
    run:
        dat = pd.read_csv(input.s, sep='\t')
        assoc = pd.read_csv(input.a, sep='\t').rename(columns={'A1':'A1_plink', 'A2':'A2_plink', 'P':'P_plink', 'OR':'OR_plink'})
        assoc.loc[:, 'SNP'] = assoc.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        adult_gwas = pd.read_csv(input.g, delim_whitespace=True).rename(columns={'SNP':'rs', 'A1':'A1_adult', 'A2':'A2_adult', 'OR':'OR_adult', 'PVAL':'PVAL_adult', 'SE':'SE_adult'})
        adult_gwas.loc[:, 'SNP'] = adult_gwas.apply(lambda row: str(row['CHR']) + ':' + str(row['BP']), axis=1)
        m1 = pd.merge(assoc, adult_gwas, on='SNP', how='left').fillna('')
        crit = m1.apply(lambda row: row['P_plink']<1e-5 or row['A1_adult']!='', axis=1)
        pd.merge(m1[crit], dat, how='left', on='SNP').to_csv(output.o, index=False, sep='\t')