# first, get variants from the parents
rule extract_variants_parents:
  input:
    aln = rules.align.output.aligned,
    idx = rules.align.output.idx,
    ref = get_reference
  output:
    var = "out/variants_parents/{sample}.tsv.gz"
  conda: "../deps/pybio/env.yml"
  container: "docker://szsctt/lr_pybio:py310"
  wildcard_constraints:
    sample = "|".join(samples.parent_name)
  shell:
    """
    python3 -m scripts.extract_features_from_sam \
     -i {input.aln} \
     -r {input.ref} \
     -o {output.var} 
    """ 

# get the location of the first and last variant in the parents
# so that we can still include reads where the alignment
# doesn't cover the full reference, as long as it covers the
# first and last variant
