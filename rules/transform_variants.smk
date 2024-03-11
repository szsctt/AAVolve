from scripts.snakemake_helpers import get_column_by_sample

# get frequency of each variant
# also get non-parental variants of high frequency. These can be 
# concatenated with parent tsv to retain these variants in subsequent steps
rule variant_frequency:
    input:
        parents = lambda wildcards: (expand("out/variants/parents/{parent}.tsv.gz", 
                                                                    parent=get_column_by_sample(wildcards, samples, "parent_name"))),
        library = "out/variants/reads/{sample}.tsv.gz"
    output:
        high_freq = "out/variants/frequency/{sample}_high.tsv.gz",
        all = "out/variants/frequency/{sample}_all.tsv.gz",
        parental = "out/variants/frequency/{sample}_parents.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    params:
        freq = lambda wildcards: f"-f {get_column_by_sample(wildcards, samples, 'non_parental_freq')}",
        input = lambda wildcards, input: f"-i {input.library}",
        parents = lambda wildcards, input: f"-p {input.parents}",
        high_freq = lambda wildcards, output: f"-o {output.high_freq}",
        all= lambda wildcards, output: f"-oa {output.all}",
        parents_out = lambda wildcards, output: f"-op {output.parental}"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.variant_frequency_long {params}
        """

# combine parental and high-frequency non-parental variants
rule combine_variants:
    input:
        parents = lambda wildcards: (expand("out/variants/parents/{parent}.tsv.gz", 
                                                                    parent=get_column_by_sample(wildcards, samples, "parent_name"))),
        high_freq = rules.variant_frequency.output.high_freq
    output:
        combined = "out/variants/combined/{sample}.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    container: "docker://szsctt/lr_pybio:py310"
    params:
        include_non_parental = lambda wildcards: get_column_by_sample(wildcards, samples, "include_non_parental"),
        combined_unzip = lambda wildcards, output: os.path.splitext(output.combined)[0]
    shell:
        """
        # cat files together but removed header from second file
        zcat {input.parents} > {params.combined_unzip}
        if [ {params.include_non_parental} == "True" ]; then
            zcat {input.high_freq} | tail -n +2 >> {params.combined_unzip}
        fi
        gzip {params.combined_unzip}
        """

# pivot long to wide to get table with one read per row
# retain all parental mutations, and discard any non-parental
rule pivot:
    input:
        library = rules.variant_frequency.input.library,
        parents = rules.combine_variants.output.combined
    output:
        pivoted_parents = "out/variants/pivot/{sample}_parents.tsv.gz",
        pivoted_seq = "out/variants/pivot/{sample}_seq.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.pivot_variants_to_wide \
         -i {input.library} \
         -p {input.parents} \
         --remove-na \
         --output-parents {output.pivoted_parents} \
         --output-seq {output.pivoted_seq}
        """

# assign parents using all columns so we have more columns for 
# unambiguous identification of parent.  
rule assign_parents:
    input:
        parents = rules.pivot.output.pivoted_parents
    output:
        assigned = "out/parents/assigned/{sample}_assigned-parents.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.assign_parents \
         -i {input.parents} \
         -o {output.assigned}
        """

# calculate frequency of each parent
rule parent_freq:
    input:
        in_file = rules.assign_parents.output.assigned
    output:
        freqs = "out/parents/freqs/{sample}_assigned-parents_freq.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.variant_frequency_wide \
            -i {input.in_file} \
            -o {output.freqs} \
            --split-counts
        """

# identify breakpoints per-read
rule ident_breakpoints:
    input:
        in_file = rules.assign_parents.output.assigned
    output:
        breakpoints = "out/parents/breaks/{sample}tsv.gz",
        break_per_read = "out/parents/breaks/{sample}-perread.tsv.gz",
        break_per_var = "out/parents/breaks/{sample}-pervar.tsv.gz"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.count_breakpoints \
            --input {input.in_file} \
            --output {output.breakpoints} \
            --summary1 {output.break_per_read} \
            --summary2 {output.break_per_var}
        """

# count distinct combinations of parents
rule distinct_parents:
    input:
        reads = rules.assign_parents.output.assigned 
    output:
        counts = "out/parents/counts/{sample}_parent-counts.tsv.gz"
    params:
        incat = lambda wildcards, input: 'zcat' if input.reads.endswith('.gz') else 'cat',
        outzip = lambda wildcards, output: 'gzip' if output.counts.endswith('.gz') else ''
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        # header
        {params.incat} {input.reads}  |\
        sed -n 1p |\
        sed 's/read_id/count/' |\
        {params.outzip} > {output.counts}

        # rest of file
        python3 -m scripts.remove_first_column -i {input.reads} |\
        sort |\
        uniq -c |\
        awk -v OFS=$'\\t' '{{$1=$1}};1' |\
        sort -t $'\\t' -k1,1nr |\
        {params.outzip} >> {output.counts}
        """

# apply variants to reference to get 'error-corrected' reads
rule apply_variants:
    input:
        variants_names_wide = rules.distinct_parents.output.counts,
        parent_variants_long = rules.combine_variants.output.combined,
        ref = get_reference
    output:
        seqs = "out/corrected/counts/{sample}_nt-seq-counts.tsv.gz"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.apply_variants \
            -v {input.variants_names_wide} \
            -p {input.parent_variants_long} \
            -r {input.ref} \
            -o {output.seqs} 
        """

rule dmat_nt:
    input:
        counts = rules.apply_variants.output.seqs
    output:
        dmat = "out/corrected/dmat/{sample}_{subset}-nt.tsv.gz",
        plot = "out/corrected/dmat/{sample}_{subset}-nt.png"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        subset = "random|first|last"
    shell:
        """
        python3 -m scripts.distance_matrix \
            --input {input.counts} \
            --output {output.dmat} \
            --plot {output.plot} \
            --distance-metric identity \
            --max-seqs 10000 \
            --selection {wildcards.subset}
        """

def get_reads_for_counting(wildcards):
    """
    Get the appropriate files for counting original reads
    If not R2C2 data, just return result of get_reads(wildcards)

    If R2C2 data, return original reads and consensus reads
    If we filtered, also return filtered reads
    """

    # for C3POa, return original fastq and processed/filtered fastq

    # if one of the parents, return parent sequences
    if wildcards.sample in parents.keys():
        return parents[wildcards.sample]
    
    # get sequencing technology
    tech = {k:v for k, v in zip(samples.name, samples.seq_technology)}[wildcards.sample]

    # if non-R2C2, this is input reads
    # if R2C2 with filtering, this is filtered reads
    # if R2C2 without filtering, this is consensus reads
    reads = [get_reads(wildcards)]

    # if not R2C2, just return reads
    if tech != 'nanopore_r2c2':
        return reads

    # for R2C2 add input reads and consensus reads
    reads.append(get_reads_for_consensus(wildcards)[0])
    reads.append("out/c3poa/{sample}/split/R2C2_Consensus.fasta.gz")
    
    return list(set(reads))
    

'''rule count_reads:
    input:
        input_reads = get_reads_for_counting, 
        variants = rules.extract_variants.output.var,
        pivoted = rules.pivot.output.pivoted_seq, 
        distinct = rules.distinct_reads.output.counts,
        distinct_nt = rules.distinct_nt.output.counts
    output:
        counts = "out/qc/{sample}_read-counts.tsv"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 scripts/count_lines.py \
         --output {output.counts} \
         --files {input}
        """
'''