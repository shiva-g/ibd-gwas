"""Plink association after imputation."""

rule plink_assoc:
    input:
        b = DATA + 'interim/bfiles_imputed_grouped/{group}/{pop}.fam',
    output:
        f = DATA + 'interim/plink_assoc/{group}/{pop}.assoc'
    singularity:
        PLINK
    log:
        LOG + 'assoc/{group}.{pop}.as'
    shell:
        "plink --bfile $(dirname {input.b})/{wildcards.pop} --assoc "
        "--out $(dirname {output.f})/{wildcards.pop} &> {log}"

rule format_assoc:
    input:
        f = DATA + 'interim/plink_assoc/{group}/{pop}.assoc'
    output:
        o = DATA + 'interim/plink_assoc_fmt/{group}/{pop}.assoc'
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
