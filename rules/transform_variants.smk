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

# assign parents using all columns  
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
        breakpoints = "out/parents/breaks/{sample}.tsv.gz",
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
rule distinct_reads:
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
        variants_names_wide = rules.distinct_reads.output.counts,
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

# translate corrected reads to amino acids
rule translate_nt:
    input:
        counts = rules.apply_variants.output.seqs
    output:
        counts = "out/corrected/counts/{sample}_aa-seq-translated.tsv.gz"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.translate_nt \
            -i {input.counts} \
            -o {output.counts}
        """

# sum counts for reads that translate to the same amino acid sequence
rule sum_nt_translated_counts:
    input:
        counts = rules.translate_nt.output.counts
    output:
        summed = "out/corrected/counts/{sample}_aa-seq-counts.tsv.gz"
    container: "docker://szsctt/lr_pybio:py310"
    params:
        cat = lambda wildcards, input: 'zcat' if input.counts.endswith('.gz') else 'cat'
    shell:
        """
        # write header
        {params.cat} {input.counts} |\
            head -n1 |\
            gzip > {output.summed}

        {params.cat} {input.counts} |\
            tail -n+2 |\
            sort -k2,2 |\
            python3 -m scripts.sum_counts |\
            sort -k1,1nr |\
            gzip >> {output.summed}
        """

def get_dmat_input(wildcards):
    if wildcards.seq_type == "nt-seq":
        return rules.apply_variants.output.seqs
    if wildcards.seq_type == "aa-seq":
        return "out/corrected/counts/{sample}_aa-seq-counts.tsv.gz"

# make distance matrix of corrected reads - either nt or aa
rule dmat:
    input:
        counts = get_dmat_input
    output:
        dmat = "out/corrected/dmat/{sample}_{subset}_{seq_type}.tsv.gz",
        plot = "out/corrected/dmat/{sample}_{subset}_{seq_type}.png"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        subset = "random|first|last",
        seq_type = "nt-seq|aa-seq"
    params:
        distance_metric = lambda wildcards: "identity" if wildcards.seq_type == "nt-seq" else "blosum62",
        max_seqs = 10000
    shell:
        """
        python3 -m scripts.distance_matrix \
            --input {input.counts} \
            --output {output.dmat} \
            --plot {output.plot} \
            --distance-metric {params.distance_metric} \
            --max-seqs {params.max_seqs} \
            --selection {wildcards.subset}
        """

# counts reads at each stage of processing
def get_reads_for_counting(wildcards):
    """
    Get the appropriate files for counting original reads
    If not R2C2 data, just return result of get_reads(wildcards)

    If R2C2 data, return original reads and consensus reads
    If we filtered, also return filtered reads
    """
    # get sequencing technology
    tech = {k:v for k, v in zip(samples.sample_name, samples.seq_tech)}[wildcards.sample]

    # if non-R2C2, this is input reads
    # if R2C2 with filtering, this is filtered reads
    # if R2C2 without filtering, this is consensus reads
    reads = [get_reads(wildcards)]

    # if not R2C2, just return reads
    if tech != 'np-cc':
        return reads

    # for R2C2 add input reads and consensus reads
    reads.append(get_column_by_sample(wildcards, samples, "read_file"))
    reads.append(rules.consensus.output.consensus_reads)
    
    return list(set(reads))
    
rule count_reads:
    input:
        input_reads = get_reads_for_counting, 
        variants = rules.extract_variants_reads.output.var,
        pivoted = rules.pivot.output.pivoted_seq, 
        distinct = rules.distinct_reads.output.counts,
        distinct_aa = rules.sum_nt_translated_counts.output.summed,
    output:
        counts = "out/qc/{sample}_read-counts.tsv"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.count_reads \
         --output {output.counts} \
         --files {input}
        """


rule report:
    input:
        counts = rules.count_reads.output.counts,
        assigned_counts = rules.distinct_reads.output.counts,
        freqs = rules.parent_freq.output.freqs,
        breaks_per_var = rules.ident_breakpoints.output.break_per_var,
        dmat_nt_first = expand(rules.dmat.output.dmat, seq_type="nt-seq", subset="first", allow_missing=True),
        dmat_aa_first = expand(rules.dmat.output.dmat, seq_type="aa-seq", subset="first", allow_missing=True),
        dmat_nt_random = expand(rules.dmat.output.dmat, seq_type="nt-seq", subset="random", allow_missing=True),
        dmat_aa_random = expand(rules.dmat.output.dmat, seq_type="aa-seq", subset="random", allow_missing=True),
    output:
        report = "out/qc/{sample}_report.html",
        tmp_notebook = "out/qc/{sample}_report.ipynb",
    container: "docker://szsctt/lr_pybio:py310"
    params:
        seq_tech = lambda wildcards: get_column_by_sample(wildcards, samples, "seq_tech"),
        tmpdir_1 = lambda wildcards, output: os.path.join(os.path.dirname(output.tmp_notebook), wildcards.sample, "runtime"),
        tmpdir_2 = lambda wildcards, output: os.path.join(os.path.dirname(output.tmp_notebook), wildcards.sample, "cache"),
        tmpdir_3 = lambda wildcards, output: os.path.join(os.path.dirname(output.tmp_notebook), wildcards.sample, "output")
    shell:
        """
        papermill scripts/report.ipynb {output.tmp_notebook} \
            -p seq_tech {params.seq_tech} \
            -p read_counts {input.counts} \
            -p assigned_parents {input.assigned_counts} \
            -p parent_frequencies {input.freqs} \
            -p breakpoints_per_var {input.breaks_per_var} \
            -p dmat_nt_first {input.dmat_nt_first} \
            -p dmat_aa_first {input.dmat_aa_first} \
            -p dmat_nt_random {input.dmat_nt_random} \
            -p dmat_aa_random {input.dmat_aa_random}

        export XDG_RUNTIME_DIR={params.tmpdir_1}
        export XDG_CACHE_HOME={params.tmpdir_2}
        export XDG_DATA_HOME={params.tmpdir_3}
        mkdir -p {params.tmpdir_1} {params.tmpdir_2} {params.tmpdir_3}
        quarto render {output.tmp_notebook} --output - > {output.report}
        """
