"""Prep hapmap for ancestry"""
# http://zzz.bwh.harvard.edu/plink/res.shtml
# https://rdrr.io/cran/plinkQC/f/inst/doc/HapMap.pdf
# http://zzz.bwh.harvard.edu/plink/dist/hapmap_CEU_r23a_filtered.zip
# http://zzz.bwh.harvard.edu/plink/dist/hapmap_YRI_r23a_filtered.zip
# http://zzz.bwh.harvard.edu/plink/dist/hapmap_JPT_CHB_r23a_filtered.zip

rule dl_hapmap:
    """hg18 hapmap phase 2"""
    input:  HTTP.remote("zzz.bwh.harvard.edu/plink/dist/hapmap_{pop}_r23a_filtered.zip", insecure=True, keep_local=True)
    output: DATA + 'raw/hapmap/{pop}.zip'
    shell:  'mv {input} {output}'

rule unzip_hapmap:
    input: DATA + 'raw/hapmap/{pop}.zip'
    output: expand(DATA + 'interim/hapmap/hapmap_{{pop}}_r23a_filtered.{s}', s=('fam', 'bim', 'bed'))
    shell:  'unzip {input} -d {DATA}/interim/hapmap/'

rule mk_coords_hapmap:
    input:
        DATA + 'interim/hapmap/hapmap_{pop}_r23a_filtered.bim'
    output:
        DATA + 'interim/hapmap/{pop}.hg18.tolift'
    shell:
        """awk '{{print "chr" $1, $4 -1, $4, $2 }}' {input} | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > {output}"""

#rule dl_chain_file:
#    input:  HTTP.remote('hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz', insecure=True, keep_local=True, allow_redirects=True)
#    output: DATA + 'raw/ucsc/hg18ToHg19.over.chain.gz'
#    shell:  'mv {input} {output}'

# rule lift_hapmap:
#     input:
#         pos = DATA + 'interim/hapmap/{pop}.tolift',
#         chain = DATA + 'raw/ucsc/hg18ToHg19.over.chain.gz'
#     output:
#         mapped = DATA + 'interim/hapmap/{pop}.mapped',
#         unmapped = DATA + 'interim/hapmap/{pop}.unmapped'
#     shell:
#         "liftOver {input.pos} {input.chain} {output.mapped} {output.unmapped}"

rule hapmap_mapped_vars:
    input:
        DATA + 'interim/hapmap/{pop}.fromhg18Tohg19.mapped',
    output:
        DATA + 'interim/hapmap/{pop}.snps'
    shell:
        "awk '{{print $4}}' {input} > {output}"

rule hapmap_mapped_pos:
    input:
        DATA + 'interim/hapmap/{pop}.mapped',
    output:
        DATA + 'interim/hapmap/{pop}.pos'
    shell:
        "awk '{{print $4, $3}}' {input} > {output}"

rule hapmap_to_hg19:
    input:
        pos = DATA + 'interim/hapmap/{pop}.pos',
        snps = DATA + 'interim/hapmap/{pop}.snps'
    output:
        expand(DATA + 'interim/hapmap/{{pop}}.{s}', s=('bed', 'bim', 'fam'))
    singularity:
        PLINK
    log:
        LOG + 'hapmap/{pop}.to_hg19'
    shell:
        "plink --bfile {DATA}interim/hapmap/hapmap_{wildcards.pop}_r23a_filtered "
        "--extract {input.snps} --update-map {input.pos} "
        "--make-bed --out {DATA}interim/hapmap/{wildcards.pop} &> {log}"

rule combine_hapmap:
    input: expand(DATA + 'interim/hapmap/{pop}.bed', pop=('CEU', 'YRI', 'JPT_CHB'))
    output:
        expand(DATA + 'interim/bfiles/hapmap.{s}', s=('bed', 'bim', 'fam'))
    singularity:
        PLINK
    log:
        LOG + 'prep/combine_hapmap'
    shell:
        "plink --merge-list {CONFIG}hapmap_bfiles "
        "--autosome --not-chr y --make-bed --out {DATA}interim/bfiles/hapmap &> {log}"
