from aavolve.snakemake_helpers import (
    get_anchor_length,
    get_column_by_sample,
    get_parents,
    get_reference_for_align,
    fill_parents,
    get_must_start_before_end_after,
)

# first, get variants from the parents
rule extract_variants_parents:
    input:
        aln = rules.align.output.aligned,
        idx = rules.align.output.idx,
        ref = lambda wildcards: get_reference_for_align(wildcards, samples)
    output:
        var = "out/variants/parents/{sample}.tsv.gz",
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.parent_name)
    params:
        anchor_length = lambda wildcards: get_anchor_length(wildcards, samples)
    log:
        "logs/extract_variants_parents/{sample}.log"
    shell:
        """
        python3 -m aavolve.extract_features_from_sam \
            -i {input.aln} \
            -r {input.ref} \
            -o {output.var} \
            --anchor-length {params.anchor_length}
        """ 

# number of parents in each parent fasta
# if there is only one we use 0 and -1 for must-start-before and must-end-after
# otherwise we use the first and last variant from rule first_last_variants_parents
rule num_parents:
    input:
        inputs_ok="out/qc/input-checks/_all.ok",
        fa = lambda wildcards: get_parents(wildcards, samples)
    output:
        num_parents = "out/variants/parents/{sample}_num_parents.txt"
    log:
        "logs/num_parents/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.parent_name)
    shell:
        """
        python3 -m aavolve.num_in_fa \
        -i {input.fa} \
        -o {output.num_parents} 
        """    

# get the location of the first and last variant in the parents
# so that we can still include reads where the alignment
# doesn't cover the full reference, as long as it covers the
# first and last variant

rule first_last_variants_parents:
    input:
        var = rules.extract_variants_parents.output.var,
        n_parents = "out/variants/parents/{sample}_num_parents.txt"
    output:
        first_last = "out/variants/parents/{sample}_first_last.txt"
    log:
        "logs/first_last_variants_parents/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.parent_name)
    shell:
        """
        NPAR=$(cat {input.n_parents})
        if [ $NPAR -eq 1 ]; then
            touch {output.first_last}
        else
            python3 -m aavolve.first_last_variant \
            -i {input.var} \
            -o {output.first_last}
        fi 
        """ 

# get variants from reads
rule extract_variants_reads:
    input:
        aln = rules.align.output.aligned,
        idx = rules.align.output.idx,
        ref = lambda wildcards: get_reference_for_align(wildcards, samples),
        first_last = lambda wildcards: fill_parents(wildcards, samples, rules.first_last_variants_parents.output.first_last),
        n_parents = lambda wildcards: fill_parents(wildcards, samples, rules.num_parents.output.num_parents)
    output:
        var = "out/variants/reads/{sample}.tsv.gz",
        read_ids = "out/variants/reads/{sample}_read-ids.txt"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    params:
        anchor_length = lambda wildcards: get_anchor_length(wildcards, samples),
        require_end_to_end_alignment = lambda wildcards: get_column_by_sample(wildcards, samples, "require_end_to_end_alignment"),
        must_start_before = lambda wildcards, input: get_must_start_before_end_after(
            first_last_file=str(input.first_last),
            n_parents_file=str(input.n_parents),
            require_end_to_end_alignment=bool(get_column_by_sample(wildcards, samples, "require_end_to_end_alignment")),
        )[0],
        must_end_after = lambda wildcards, input: get_must_start_before_end_after(
            first_last_file=str(input.first_last),
            n_parents_file=str(input.n_parents),
            require_end_to_end_alignment=bool(get_column_by_sample(wildcards, samples, "require_end_to_end_alignment")),
        )[1],
    log:
        "logs/extract_variants_reads/{sample}.log"
    shell:
        """
        echo "Using must_start_before={params.must_start_before} and must_end_after={params.must_end_after} (require_end_to_end_alignment={params.require_end_to_end_alignment})."

        python3 -m aavolve.extract_features_from_sam \
            -i {input.aln} \
            -r {input.ref} \
            -o {output.var} \
            -O {output.read_ids} \
            --must-start-before {params.must_start_before} \
            --must-end-after {params.must_end_after} \
            --anchor-length {params.anchor_length}
        """

# count the number of reads we ended up with data for (some will be excluded
# because they didn't cover enough of the reference)
rule count_variant_reads:
    input:
        var = rules.extract_variants_reads.output.var
    output:
        read_count = "out/variants/reads/{sample}_read-count.txt"
    params:
        cat = lambda wildcards, input: "zcat" if input.var.endswith('.gz') else 'cat'
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    log:
        "logs/count_variant_reads/{sample}.log"
    shell:
        """
        {params.cat} {input.var} | cut -f3 -d$'\\t' | uniq | wc -l > {output.read_count}
        """
