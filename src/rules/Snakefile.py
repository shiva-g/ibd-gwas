"""Main snakefile"""

include: "const.py"
include: "sf_prep.py"

rule all:
    output: LOG + 'DONE'
    shell:  'touch {output}'
