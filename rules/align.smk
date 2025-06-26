import numpy as np
from aavolve.snakemake_helpers import get_reads, get_reference, minimap2_params_with_default

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
    params:
        minimap2_params = lambda wildcards: minimap2_params_with_default(wildards, samples),
    shell:
        """
        minimap2 -t {threads} -a {params.minimap2_params} {input.reference} {input.reads} --MD |\
            samtools sort -o {output.aligned} -
        
        samtools index {output.aligned}
        """