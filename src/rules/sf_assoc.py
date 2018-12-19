"""Plink association after imputation."""

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
        df = df[pd.notnull(df['P'])]
        df.to_csv(output.o, index=False, sep='\t')

# http://www.gettinggeneticsdone.com/2011/04/annotated-manhattan-plots-and-qq-plots.html
rule assoc_plot:
    input:
        DATA + 'interim/plink_assoc_fmt/{group}/eur.assoc'
    output:
        PLOTS + 'manhattan.{group}.png',
        PLOTS + 'qq.{group}.png'
    singularity:
        "docker://manninglab/metal"
    shell:
        "Rscript {SCRIPTS}plot_manhattan.R {input} {output}"
