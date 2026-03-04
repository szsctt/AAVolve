import os

from aavolve.snakemake_helpers import (
    get_column_by_sample,
    get_reference_for_align,
    is_fastq,
    get_reads_for_counting,
    format_input_reads,
    get_dmat_input,
)
from aavolve.snakemake_helpers import get_non_parental_high_freq_variants_for_sample
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
    log:
        "logs/variant_frequency/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.variant_frequency_long {params}
        """

rule combine_non_parental_variants_group:
    """
    Combine high-frequency non-parental variants across samples that share the same
    (parent_file, reference_file) pair and have include_non_parental=True, and (if configured)
    share the same non_parental_group.
    """
    input:
        variants=lambda wildcards: expand(
            "out/variants/frequency/{sample}_high.tsv.gz",
            sample=non_parental_variant_groups[wildcards.group_id],
        ),
    output:
        combined="out/variants/frequency/groups/{group_id}_high.tsv.gz",
    log:
        "logs/combine_non_parental_variants_group/{group_id}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.combine_non_parental_variants \
            --output {output.combined} \
            {input.variants} \
            2> {log}
        """

# combine parental and high-frequency non-parental variants
rule combine_variants:
    input:
        parents = lambda wildcards: (expand("out/variants/parents/{parent}.tsv.gz", 
                                                                    parent=get_column_by_sample(wildcards, samples, "parent_name"))),
        high_freq = lambda wildcards: get_non_parental_high_freq_variants_for_sample(
            wildcards,
            samples,
            sample_to_non_parental_group_id,
            non_parental_variant_groups,
        )
    output:
        combined = "out/variants/combined/{sample}.tsv.gz"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    container: "docker://szsctt/lr_pybio:py310"
    params:
        include_non_parental = lambda wildcards: get_column_by_sample(wildcards, samples, "include_non_parental"),
        combined_unzip = lambda wildcards, output: os.path.splitext(output.combined)[0]
    log:
        "logs/combine_variants/{sample}.log"
    shell:
        """
        # cat files together but removed header from second file
        zcat {input.parents} > {params.combined_unzip}
        if [ {params.include_non_parental} == "True" ]; then
            zcat {input.high_freq} | tail -n +2 >> {params.combined_unzip}
        fi
        gzip {params.combined_unzip}
        """

rule warn_non_parental_outside_window:
    input:
        high_freq = rules.variant_frequency.output.high_freq,
        first_last = lambda wildcards: expand(
            "out/variants/parents/{parent}_first_last.txt",
            parent=get_column_by_sample(wildcards, samples, "parent_name"),
        ),
        n_parents = lambda wildcards: expand(
            "out/variants/parents/{parent}_num_parents.txt",
            parent=get_column_by_sample(wildcards, samples, "parent_name"),
        ),
    output:
        warning = "out/qc/warnings/{sample}_non_parental_outside_window.txt"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    params:
        include_non_parental = lambda wildcards: get_column_by_sample(wildcards, samples, "include_non_parental"),
        require_end_to_end_alignment = lambda wildcards: get_column_by_sample(wildcards, samples, "require_end_to_end_alignment"),
    log:
        "logs/warn_non_parental_outside_window/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        NPAR=$(cat {input.n_parents})
        python3 -m aavolve.warn_non_parental_outside_window \
            --high-freq-variants {input.high_freq} \
            --parent-first-last {input.first_last} \
            --n-parents $NPAR \
            --include-non-parental {params.include_non_parental} \
            --require-end-to-end-alignment {params.require_end_to_end_alignment} \
            --output {output.warning} \
            2> {log}
        """


rule warn_parents_dropped:
    input:
        bam=lambda wildcards: f"out/aligned/{get_column_by_sample(wildcards, samples, 'parent_name')}.bam",
        bai=lambda wildcards: f"out/aligned/{get_column_by_sample(wildcards, samples, 'parent_name')}.bam.bai",
        parents=lambda wildcards: {k: v for k, v in zip(samples.parent_name, samples.parent_file)}[
            get_column_by_sample(wildcards, samples, "parent_name")
        ],
        reference=lambda wildcards: get_reference_for_align(
            type(
                "WC",
                (),
                {"sample": get_column_by_sample(wildcards, samples, "parent_name")},
            )(),
            samples,
        ),
    output:
        warning="out/qc/warnings/{sample}_parents_dropped.txt",
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    log:
        "logs/warn_parents_dropped/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.warn_parents_dropped \
            --parents-fasta {input.parents} \
            --bam {input.bam} \
            --reference-fasta {input.reference} \
            --output {output.warning} \
            2> {log}
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
        max_group_distance = lambda wildcards: '--max-distance-frac 0' if not get_column_by_sample(wildcards, samples, "group_vars") else f'--max-distance-frac {get_column_by_sample(wildcards, samples, "max_group_distance")}',
    log:
        "logs/pivot/{sample}.log"
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
    log:
        "logs/assign_parents/{sample}.log"
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
    log:
        "logs/parent_freq/{sample}.log"
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
    log:
        "logs/ident_breakpoints/{sample}.log"
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
        counts = "out/parents/counts/{sample}_parent-counts.tsv.gz",
        members = "out/parents/counts/{sample}_parent-counts-members.tsv.gz"
    log:
        "logs/distinct_reads/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.distinct_reads \
            -i {input.reads} \
            -o {output.counts} \
            -m {output.members}
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
    params:
        group_variants = lambda wildcards: '--group-vars' if get_column_by_sample(wildcards, samples, "group_vars") else '',
        group_dist = lambda wildcards: f'--group-dist {get_column_by_sample(wildcards, samples, "group_vars_dist")}' if get_column_by_sample(wildcards, samples, "group_vars_dist") else '',
    log:
        "logs/apply_variants/{sample}.log"
    shell:
        """
        python3 -m aavolve.apply_variants \
            -v {input.variants_names_wide} \
            -p {input.parent_variants_long} \
            -r {input.ref} \
            -o {output.seqs}  \
            {params}
        """

# translate corrected reads to amino acids
rule translate_nt:
    input:
        counts = rules.apply_variants.output.seqs
    output:
        counts = temp("out/corrected/counts/{sample}_aa-seq-translated.tsv.gz")
    log:
        "logs/translate_nt/{sample}.log"
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
    log:
        "logs/sum_nt_translated_counts/{sample}.log"
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
        max_seqs = config.get("max_seqs", MAX_SEQS)
    log:
        "logs/dmat/{sample}_{subset}_{seq_type}.log"
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
    log:
        "logs/count_reads/{sample}.log"
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
        variant_freq_all = rules.variant_frequency.output.all,
        combined_variants = rules.combine_variants.output.combined,
        breaks_per_var = rules.ident_breakpoints.output.break_per_var,
        variant_window_warning = rules.warn_non_parental_outside_window.output.warning,
        parents_dropped_warning = rules.warn_parents_dropped.output.warning,
        coverage_depth = rules.coverage_depth.output.depth,
        pretrim_msa = "out/qc/trimming/{sample}.pretrim.mafft.fasta",
        posttrim_msa = lambda wildcards: (
            f"out/qc/trimming/{wildcards.sample}.mafft.fasta"
            if ("trim" in samples.columns and get_column_by_sample(wildcards, samples, "trim"))
            else []
        ),
        dmat_nt_first = expand(rules.dmat.output.dmat, seq_type="nt-seq", subset="first", allow_missing=True),
        dmat_aa_first = expand(rules.dmat.output.dmat, seq_type="aa-seq", subset="first", allow_missing=True),
        dmat_nt_random = expand(rules.dmat.output.dmat, seq_type="nt-seq", subset="random", allow_missing=True),
        dmat_aa_random = expand(rules.dmat.output.dmat, seq_type="aa-seq", subset="random", allow_missing=True),
        report_template = lambda wildcards: os.path.join(workflow.basedir, "aavolve/report.ipynb")
    output:
        report = "out/reports/sample_reports/{sample}_report.html",
        tmp_notebook = temp("out/reports/sample_reports/{sample}_report.ipynb"),
    log:
        "logs/report/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    params:
        seq_tech = lambda wildcards: get_column_by_sample(wildcards, samples, "seq_tech"),
        trim = lambda wildcards: get_column_by_sample(wildcards, samples, "trim") if "trim" in samples.columns else False,
        adapter_5 = lambda wildcards: (
            "" if ("adapter_5" not in samples.columns or str(get_column_by_sample(wildcards, samples, "adapter_5")).lower() == "nan")
            else str(get_column_by_sample(wildcards, samples, "adapter_5"))
        ),
        adapter_3 = lambda wildcards: (
            "" if ("adapter_3" not in samples.columns or str(get_column_by_sample(wildcards, samples, "adapter_3")).lower() == "nan")
            else str(get_column_by_sample(wildcards, samples, "adapter_3"))
        ),
        anchors = lambda wildcards: int(get_column_by_sample(wildcards, samples, "anchors")) if "anchors" in samples.columns else 0,
        include_non_parental = lambda wildcards: get_column_by_sample(wildcards, samples, "include_non_parental"),
        non_parental_freq = lambda wildcards: get_column_by_sample(wildcards, samples, "non_parental_freq") if "non_parental_freq" in samples.columns else 0.0,
        require_end_to_end_alignment = lambda wildcards: get_column_by_sample(wildcards, samples, "require_end_to_end_alignment") if "require_end_to_end_alignment" in samples.columns else False,
        group_vars = lambda wildcards: get_column_by_sample(wildcards, samples, "group_vars") if "group_vars" in samples.columns else "",
        group_vars_dist = lambda wildcards: get_column_by_sample(wildcards, samples, "group_vars_dist") if "group_vars_dist" in samples.columns else "",
        max_group_distance = lambda wildcards: get_column_by_sample(wildcards, samples, "max_group_distance") if "max_group_distance" in samples.columns else "",
        minimap2_params = lambda wildcards: get_column_by_sample(wildcards, samples, "minimap2_params") if "minimap2_params" in samples.columns else "",
        report_basename = lambda wildcards, output: os.path.basename(output.tmp_notebook),
        report_dir = lambda wildcards, output: os.path.dirname(output.report),
        deno_mem_mb = lambda wildcards: int(config.get("quarto_deno_max_old_space_size_mb", 4096)),
    shell:
        """
        pwd
        papermill {input.report_template} {output.tmp_notebook} \
            -p seq_tech {params.seq_tech} \
            -p read_counts {input.counts} \
            -p assigned_parents {input.assigned_counts} \
            -p parent_frequencies {input.freqs} \
            -p variant_freq_all {input.variant_freq_all} \
            -p combined_variants {input.combined_variants} \
            -p breakpoints_per_var {input.breaks_per_var} \
            -p variant_window_warning {input.variant_window_warning} \
            -p parents_dropped_warning {input.parents_dropped_warning} \
            -p coverage_depth {input.coverage_depth} \
            -p pretrim_msa {input.pretrim_msa} \
            -p posttrim_msa "{input.posttrim_msa}" \
            -p dmat_nt_first {input.dmat_nt_first} \
            -p dmat_aa_first {input.dmat_aa_first} \
            -p dmat_nt_random {input.dmat_nt_random} \
            -p dmat_aa_random {input.dmat_aa_random} \
            -p trim {params.trim} \
            -p adapter_5 "{params.adapter_5}" \
            -p adapter_3 "{params.adapter_3}" \
            -p anchors {params.anchors} \
            -p include_non_parental {params.include_non_parental} \
            -p non_parental_freq {params.non_parental_freq} \
            -p require_end_to_end_alignment {params.require_end_to_end_alignment} \
            -p group_vars {params.group_vars} \
            -p group_vars_dist {params.group_vars_dist} \
            -p max_group_distance {params.max_group_distance} \
            -p minimap2_params "{params.minimap2_params}"

        cd {params.report_dir}
        unset QUARTO_DENO DENO
        export QUARTO_DENO_EXTRA_OPTIONS="--v8-flags=--max-old-space-size={params.deno_mem_mb}"
        quarto render {params.report_basename}
        """
