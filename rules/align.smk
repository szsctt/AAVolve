from aavolve.snakemake_helpers import (
    get_anchor_length,
    get_anchor_reads_suffix,
    get_anchor_seed,
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
    input:
        inputs_ok="out/qc/input-checks/_all.ok"
    output:
        anchors="out/anchors/{sample}.fasta"
    params:
        length=lambda wildcards: get_anchor_length(wildcards, samples),
        seed=lambda wildcards: get_anchor_seed(wildcards, samples)
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
        inputs_ok="out/qc/input-checks/_all.ok",
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
        inputs_ok="out/qc/input-checks/_all.ok",
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

rule trim_reads:
    """
    Trim adapters from reads using cutadapt when configured (trim=True).
    Uses linked-adapter syntax '<5>...<3>' and retains only trimmed reads (--discard-untrimmed).
    """
    input:
        inputs_ok="out/qc/input-checks/_all.ok",
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

rule qc_trimming_msa:
    """
    QC adapter trimming by aligning the reference plus the first 200 trimmed reads.
    """
    input:
        inputs_ok="out/qc/input-checks/_all.ok",
        trimmed="out/trimmed/{sample}.trimmed.gz",
        reference=lambda wildcards: get_reference(wildcards, samples),
    output:
        msa="out/qc/trimming/{sample}.mafft.fasta",
    threads: 4
    log:
        "logs/qc_trimming_msa/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.qc_trimming_msa \
            --reference {input.reference} \
            --trimmed {input.trimmed} \
            --output {output.msa} \
            --max-seqs 200 \
            --threads {threads} \
            2> {log}
        """

rule qc_pretrim_msa:
    """
    QC whether reads need trimming by aligning the reference plus the first 200 untrimmed reads.
    """
    input:
        inputs_ok="out/qc/input-checks/_all.ok",
        reads=lambda wildcards: get_reads(wildcards, samples),
        reference=lambda wildcards: get_reference(wildcards, samples),
    output:
        msa="out/qc/trimming/{sample}.pretrim.mafft.fasta",
    threads: 4
    log:
        "logs/qc_pretrim_msa/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.qc_trimming_msa \
            --reference {input.reference} \
            --reads {input.reads} \
            --output {output.msa} \
            --max-seqs 200 \
            --threads {threads} \
            2> {log}
        """


# map to one of the parental references.  The choice of reference is arbitrary
rule align:
    input:
        inputs_ok="out/qc/input-checks/_all.ok",
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
