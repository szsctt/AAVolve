from aavolve.snakemake_helpers import get_column_by_sample, is_fastq, get_reads_for_counting, format_input_reads, get_dmat_input
from aavolve.utils import MAX_SEQS

# get frequency of each variant
# also get non-parental variants of high frequency. These can be 
# concatenated with parent tsv to retain these variants in subsequent steps
rule variant_frequency:
    input:
        parents = lambda wildcards: (expand("out/variants/parents/{parent}.tsv.gz", 
                                                                    parent=get_column_by_sample(wildcards, samples, "parent_name"))),
        library = rules.extract_variants_reads.output.var,
        library_read_ids = rules.extract_variants_reads.output.read_ids,
    output:
        high_freq = "out/variants/frequency/{sample}_high.tsv.gz",
        all = "out/variants/frequency/{sample}_all.tsv.gz",
        parental = "out/variants/frequency/{sample}_parents.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name),
    params:
        freq = lambda wildcards: f"-f {get_column_by_sample(wildcards, samples, 'non_parental_freq')}",
        input = lambda wildcards, input: f"-i {input.library}",
        input_read_ids = lambda wildcards, input: f"-r {input.library_read_ids}",
        parents = lambda wildcards, input: f"-p {input.parents}",
        high_freq = lambda wildcards, output: f"-o {output.high_freq}",
        all= lambda wildcards, output: f"-oa {output.all}",
        parents_out = lambda wildcards, output: f"-op {output.parental}"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.variant_frequency_long {params}
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
rule pivot:
    input:
        library = rules.variant_frequency.input.library,
        library_read_ids = rules.extract_variants_reads.output.read_ids,
        parents = rules.combine_variants.output.combined
    output:
        pivoted_parents = "out/variants/pivot/{sample}_parents.tsv.gz",
        pivoted_seq = "out/variants/pivot/{sample}_seq.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    params:
        group_variants = lambda wildcards: '--group-vars' if get_column_by_sample(wildcards, samples, "group_vars") else '',
        group_dist = lambda wildcards: f'--group-dist {get_column_by_sample(wildcards, samples, "group_vars_dist")}' if get_column_by_sample(wildcards, samples, "group_vars_dist") else '',
        max_group_distance = lambda wildcards: '--max-group-distance 0' if not get_column_by_sample(wildcards, samples, "group_vars") else f'--max-group-distance {get_column_by_sample(wildcards, samples, "max_group_distance")}',

    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.pivot_variants_to_wide \
         -i {input.library} \
         -r {input.library_read_ids} \
         -p {input.parents} \
         --remove-na \
         --output-parents {output.pivoted_parents} \
         --output-seq {output.pivoted_seq} \
         {params}
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
        python3 -m aavolve.assign_parents \
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
        python3 -m aavolve.variant_frequency_wide \
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
        python3 -m aavolve.count_breakpoints \
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
        python3 -m aavolve.remove_first_column -i {input.reads} |\
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
        ref = lambda wildcards: get_reference(wildcards, samples)
    output:
        seqs = "out/corrected/counts/{sample}_nt-seq-counts.tsv.gz"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.apply_variants \
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
        python3 -m aavolve.translate_nt \
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
            awk 'NR==1' |\
            gzip > {output.summed}

        {params.cat} {input.counts} |\
            awk 'NR!=1' |\
            sort -k2,2 |\
            python3 -m aavolve.sum_counts |\
            sort -k1,1nr |\
            gzip >> {output.summed}
        """

# make distance matrix of corrected reads - either nt or aa
rule dmat:
    input:
        counts = lambda wildcards: get_dmat_input(wildcards, rules.apply_variants.output.seqs, "out/corrected/counts/{sample}_aa-seq-counts.tsv.gz")
    output:
        dmat = "out/corrected/dmat/{sample}_{subset}_{seq_type}.tsv.gz",
        plot = "out/corrected/dmat/{sample}_{subset}_{seq_type}.png"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        subset = "random|first|last",
        seq_type = "nt-seq|aa-seq"
    params:
        distance_metric = lambda wildcards: "identity" if wildcards.seq_type == "nt-seq" else "blosum62",
        max_seqs = MAX_SEQS
    shell:
        """
        python3 -m aavolve.distance_matrix \
            --input {input.counts} \
            --output {output.dmat} \
            --plot {output.plot} \
            --distance-metric {params.distance_metric} \
            --max-seqs {params.max_seqs} \
            --selection {wildcards.subset}
        """

    
rule count_reads:
    input:
        input_reads = lambda wildcards: get_reads_for_counting(wildcards, samples, rules.consensus.output.consensus_reads, rules.filter_consensus.output.filt), 
        variants = rules.extract_variants_reads.output.read_ids,
        pivoted = rules.pivot.output.pivoted_seq, 
        distinct = rules.distinct_reads.output.counts,
        distinct_aa = rules.sum_nt_translated_counts.output.summed,
    output:
        counts = "out/qc/{sample}_read-counts.tsv"
    container: "docker://szsctt/lr_pybio:py310"
    params:
        input_file = lambda wildcards, input: format_input_reads(input.input_reads),
        variants = lambda wildcards, input: f"--variant-read-ids {input.variants}",
        pivoted = lambda wildcards, input: f"--pivoted-tsv-files {input.pivoted}",
        distinct = lambda wildcards, input: f"--distinct-read-counts-files {input.distinct} {input.distinct_aa}",
    shell:
        """
        python3 -m aavolve.count_reads \
         --output {output.counts} \
         {params}
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
        report_basename = lambda wildcards, output: os.path.basename(output.tmp_notebook)
    shell:
        """
        papermill aavolve/report.ipynb {output.tmp_notebook} \
            -p seq_tech {params.seq_tech} \
            -p read_counts {input.counts} \
            -p assigned_parents {input.assigned_counts} \
            -p parent_frequencies {input.freqs} \
            -p breakpoints_per_var {input.breaks_per_var} \
            -p dmat_nt_first {input.dmat_nt_first} \
            -p dmat_aa_first {input.dmat_aa_first} \
            -p dmat_nt_random {input.dmat_nt_random} \
            -p dmat_aa_random {input.dmat_aa_random}

        cd out/qc
        quarto render {params.report_basename}
        """

'''

        #export XDG_RUNTIME_DIR={params.tmpdir_1}
        #export XDG_CACHE_HOME={params.tmpdir_2}
        #export XDG_DATA_HOME={params.tmpdir_3}
        #mkdir -p {params.tmpdir_1} {params.tmpdir_2} {params.tmpdir_3}
        #quarto render {output.tmp_notebook} --output - > {output.report}'''