"""Split samples by SNP"""

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
        DATA + 'interim/extracted_snps/{pop}.ped'
    singularity:
        PLINK
    log:
        LOG + 'extracted_snps/{pop}'
    shell:
        'plink --bfile $(dirname {input.f})/{wildcards.pop} '
        '--extract {input.s} --recode '
        '--out $(dirname {output})/{wildcards.pop} &> {log}'
