from aavolve.snakemake_helpers import (
    get_anchor_length,
    get_anchor_reads_suffix,
    get_linked_adapters_for_sample,
    get_reads,
    get_reads_for_align,
    get_reads_for_anchor_input,
    get_reference,
    get_reference_for_align,
    get_anchored_reference_path,
    minimap2_params_with_default,
)

rule generate_anchor_sequences:
    output:
        anchors="out/anchors/{sample}.fasta"
    params:
        length=lambda wildcards: get_anchor_length(wildcards, samples),
        seed=lambda wildcards: wildcards.sample
    log:
        "logs/generate_anchor_sequences/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.generate_anchors \
            --length {params.length} \
            --output {output.anchors} \
            --seed "{params.seed}"
        """


rule anchor_reads:
    input:
        anchors="out/anchors/{sample}.fasta",
        reads=lambda wildcards: get_reads_for_anchor_input(wildcards, samples)
    output:
        anchored="out/anchors/reads/{sample}{anchor_suffix}"
    params:
        expected_suffix=lambda wildcards: get_anchor_reads_suffix(wildcards, samples)
    log:
        "logs/anchor_reads/{sample}{anchor_suffix}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.apply_anchors \
            --input {input.reads} \
            --output {output.anchored} \
            --anchors {input.anchors}
        """


rule anchor_reference:
    input:
        anchors="out/anchors/{sample}.fasta",
        reference=lambda wildcards: get_reference(wildcards, samples)
    output:
        anchored="out/anchors/references/{sample}.fasta"
    log:
        "logs/anchor_reference/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.apply_anchors \
            --input {input.reference} \
            --output {output.anchored} \
            --anchors {input.anchors}
        """

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
    log:
        "logs/trim_reads/{sample}.log"
    container: "docker://quay.io/biocontainers/cutadapt:5.1--py39hbcbf7aa_0"
    shell:
        """
        # linked adapters validated in params function; if none, this rule should not be required
        (
            cutadapt -j {threads} -a '{params.linked}' --discard-untrimmed -o - {input.reads} \
                | pigz -p {threads} > {output.trimmed}
        ) 2> {log}
        """


rule align:
    input:
        reads = lambda wildcards: get_reads_for_align(wildcards, samples),
        reference = lambda wildcards: get_reference_for_align(wildcards, samples)
    output:
        aligned = "out/aligned/{sample}.bam",
        idx = "out/aligned/{sample}.bam.bai",
    container: "docker://szsctt/lr_align"
    threads: 8
    params:
        minimap2_params = lambda wildcards: minimap2_params_with_default(wildcards, samples),
    log:
        "logs/align/{sample}.log"
    shell:
        """
        minimap2 -t {threads} -a {params.minimap2_params} {input.reference} {input.reads} --MD |\
            samtools sort -o {output.aligned} -
        
        samtools index {output.aligned}
        """