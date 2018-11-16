"""Main snakefile"""

include: "const.py"
include: "sf_prep.py"
include: "sf_hapmap.py"

rule all:
    output: LOG + 'DONE'
    shell:  'touch {output}'
