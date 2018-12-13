"""Run plink association."""

rule plink_assoc:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/eur.fam',
    output:
        f = DATA + 'interim/plink_assoc/{group}/eur.assoc'
    singularity:
        PLINK
    log:
        LOG + 'assoc/{group}.as'
    shell:
        "plink --bfile $(dirname {input.b})/eur --assoc "
        "--out $(dirname {output.f})/eur &> {log}"

rule format_assoc:
    input:
        f = DATA + 'interim/plink_assoc/{group}/eur.assoc'
    output:
        o = DATA + 'interim/plink_assoc_fmt/{group}/eur.assoc'
    run:
        df = pd.read_csv(input.f, delim_whitespace=True).sort_values(by='CHISQ', ascending=False)
        df.to_csv(output.o, index=False, sep='\t')
