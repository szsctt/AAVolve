from aavolve.snakemake_helpers import (
    get_linked_adapters_for_sample,
    get_reads,
    get_reads_for_align,
    get_reference,
    minimap2_params_with_default,
)

# map to one of the parental references.  The choice of reference is arbitrary
rule trim_reads:
    """
    Trim adapters from reads using cutadapt when configured (trim=True).
    Uses linked-adapter syntax '<5>...<3>' and retains only trimmed reads (--discard-untrimmed).
    """
    input:
        reads = lambda wildcards: get_reads(wildcards, samples)
    output:
        trimmed = "out/trimmed/{sample}.trimmed.gz"
    params:
        linked = lambda wildcards: get_linked_adapters_for_sample(wildcards, samples)
    threads: 4
    container: "docker://quay.io/biocontainers/cutadapt:5.1--py39hbcbf7aa_0"
    shell:
        """
        # linked adapters validated in params function; if none, this rule should not be required
        cutadapt -j {threads} -a '{params.linked}' --discard-untrimmed -o - {input.reads} | pigz -p {threads} > {output.trimmed}
        """


rule align:
    input:
        reads = lambda wildcards: get_reads_for_align(wildcards, samples),
        reference = lambda wildcards: get_reference(wildcards, samples)
    output:
        aligned = "out/aligned/{sample}.bam",
        idx = "out/aligned/{sample}.bam.bai",
    conda: "../deps/align/env.yml"
    container: "docker://szsctt/lr_align"
    threads: 8
    params:
        minimap2_params = lambda wildcards: minimap2_params_with_default(wildcards, samples),
    shell:
        """
        minimap2 -t {threads} -a {params.minimap2_params} {input.reference} {input.reads} --MD |\
            samtools sort -o {output.aligned} -
        
        samtools index {output.aligned}
        """