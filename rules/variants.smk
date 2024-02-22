from scripts.snakemake_helpers import get_column_by_sample

# first, get variants from the parents
rule extract_variants_parents:
  input:
    aln = rules.align.output.aligned,
    idx = rules.align.output.idx,
    ref = get_reference
  output:
    var = "out/variants/parents/{sample}.tsv.gz"
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

def get_parents(wildcards):
    parents = {k:v for k,v in zip(samples.parent_name, samples.parent_file)}
    return parents[wildcards.sample]

# number of parents in each parent fasta
# if there is only one we use 0 and -1 for must-start-before and must-end-after
# otherwise we use the first and last variant from rule first_last_variants_parents
rule num_parents:
    input:
        fa = get_parents
    output:
        num_parents = "out/variants/parents/{sample}_num_parents.txt"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.parent_name)
    shell:
        """
        python3 -m scripts.num_in_fa \
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
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.parent_name)
    shell:
        """
        NPAR=$(cat {input.n_parents})
        if [ $NPAR -eq 1 ]; then
            touch {output.first_last}
        else
            python3 -m scripts.first_last_variant \
            -i {input.var} \
            -o {output.first_last}
        fi 
        """ 


def get_parents_first_last(wildcards):
    parent = get_column_by_sample(wildcards, samples, "parent_name")
    return f"out/variants/parents/{parent}_first_last.txt"

def get_num_parents(wildcards):
    parent = get_column_by_sample(wildcards, samples, "parent_name")
    return f"out/variants/parents/{parent}_num_parents.txt"

# get variants from reads
rule extract_variants_reads:
    input:
        aln = rules.align.output.aligned,
        idx = rules.align.output.idx,
        ref = get_reference,
        first_last = get_parents_first_last,
        n_parents = get_num_parents
    output:
        var = "out/variants/reads/{sample}.tsv.gz"
    conda: "../deps/pybio/env.yml"
    container: "docker://szsctt/lr_pybio:py310"
    wildcard_constraints:
        sample = "|".join(samples.sample_name)
    shell:
        """
        NPAR=$(cat {input.n_parents})
        if [ $NPAR -eq 1 ]; then
            FIRST=0
            LAST=-1
        else
            FIRST=$(sed -n 1p {input.first_last})
            LAST=$(sed -n 2p {input.first_last})
        fi
        echo "For $NPAR parents, using first $FIRST and last $LAST as start and end."

        python3 -m scripts.extract_features_from_sam \
            -i {input.aln} \
            -r {input.ref} \
            -o {output.var} \
            --must-start-before $FIRST \
            --must-end-after $LAST \
        """