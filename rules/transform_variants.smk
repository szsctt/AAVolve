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
        parents = rules.variant_frequency.input.parents
    output:
        pivoted_parents = "out/variants/pivot/{sample}_parents.tsv.gz",
        pivoted_seq = "out/variants/pivot/{sample}_seq.tsv.gz"
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
# unambiguous identification of parent.  We ignore NA in this case,
# assuming them to be errors.  Both reads with NA and non-aa changing cols
# are removed in the next step
rule assign_parents:
    input:
        parents = rules.pivot.output.pivoted_parents
    output:
        assigned = "out/parents/{sample}_assigned-parents.tsv.gz"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.assign_parents \
         -i {input.parents} \
         -o {output.assigned}
        """

