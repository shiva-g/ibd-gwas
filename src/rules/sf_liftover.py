from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule dl_chain_file:
    input:
        HTTP.remote('hgdownload.cse.ucsc.edu/goldenPath/{from}/liftOver/{from}To{to}.over.chain.gz', insecure=True, keep_local=True, allow_redirects=True)
    output:
        '{chain_prefix}ucsc/chains/{from}To{to}.over.chain.gz'
    shell:
        'mv {input} {output}'

# rule lift_hapmap:
#    input:
#        pos = '{bed_prefix}.{from}.tolift',
#        chain = '{chain_prefix}ucsc/chains/{from}To{to}.over.chain.gz'
#    output:
#        mapped = '{liftoot}.from{from}To{to}.mapped',
#        unmapped = '{liftout}.from{from}To{to}.unmapped'
#    shell:
#        "liftOver {input.pos} {input.chain} {output.mapped} {output.unmapped}"
