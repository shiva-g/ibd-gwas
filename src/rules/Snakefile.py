"""Main snakefile"""

include: "const.py"
include: "sf_manifest.py"
include: "sf_filter.py"
include: "sf_prep.py"
include: "sf_hapmap.py"
include: "sf_mds.py"

rule all:
    output: LOG + 'DONE'
    shell:  'touch {output}'
