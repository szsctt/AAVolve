# get frequency of each variant
# also get non-parental variants of high frequency. These can be 
# concatenated with parent tsv to retain these variants in subsequent steps
rule variant_frequency:
  input:
    parents = lambda wildcards: expand("out/variants/{parent}.tsv.gz", parent=get_parents(wildcards)),
    library = "out/variants/{sample}.tsv.gz"
  output:
    high_freq = "out/variants/frequency/{sample}.tsv.gz",
    all = "out/variants/frequency/{sample}_all.tsv.gz",
    parental = "out/variants/frequency/{sample}_parents.tsv.gz"
  params:
    freq = "-f 0.1",
    input = lambda wildcards, input: f"-i {input.library}",
    parents = lambda wildcards, input: f"-p {input.parents}",
    high_freq = lambda wildcards, output: f"-o {output.high_freq}",
    all= lambda wildcards, output: f"-oa {output.all}",
    parents_out = lambda wildcards, output: f"-op {output.parental}"
  conda: "envs/py37.yml"
  container: "docker://szsctt/lr_pybio:py310"
  shell:
    """
    python3 scripts/variant_frequency_long.py {params}
    """