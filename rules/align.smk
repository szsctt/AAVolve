import numpy as np
from aavolve.snakemake_helpers import get_reads, get_reference

# map to one of the parental references.  The choice of reference is arbitrary
rule align:
    input:
        reads = lambda wildcards: get_reads(wildcards, samples),
        reference = lambda wildcards: get_reference(wildcards, samples)
    output:
        aligned = "out/aligned/{sample}.bam",
        idx = "out/aligned/{sample}.bam.bai",
    conda: "../deps/align/env.yml"
    container: "docker://szsctt/lr_align"
    threads: 8
    shell:
        """
        minimap2 -t {threads} -ax map-hifi {input.reference} {input.reads} -B 1.5 --end-bonus 5 --MD |\
            samtools sort -o {output.aligned} -
        
        samtools index {output.aligned}
        """