from scripts.snakemake_helpers import get_column_by_sample

def get_reads(wildcards):
    """
    Get appropriate reads for wildcards.sample
    Either parental references,
    output from C3POa for nanopore R2C2 reads
    or just fastq otherwise
    """
    # make a dictionary of parents
    parents = {}
    for k, v in zip(samples['parent_name'], samples['parent_file']):
        parents[k] = v
    
    # if one of the parents, return parent sequences
    if wildcards.sample in parents.keys():
        return parents[wildcards.sample]
    
    # get sequencing technology
    tech = get_column_by_sample(wildcards, samples, 'seq_tech')
   
    # if nanopore r2c2, return consensus reads
    if tech == 'nanopore_r2c2':
        
        # check if we want to filter for repeats
        repeats = get_column_by_sample(wildcards, samples, 'min_reps')

        if math.isnan(repeats):
            return f"out/c3poa/{wildcards.sample}/split/R2C2_Consensus.fasta.gz"

        else:
            return f"out/c3poa_filt/{wildcards.sample}.fasta.gz"

    # otherwise, just return reads
    return get_column_by_sample(wildcards, samples, 'read_file')

def get_reference(wildcards):
    """
    Get appropriate reference for wildcards.sample
    Either parental references,
    or just the reference otherwise
    """
    # make a dictionary of parents
    parents = {}
    for k, v in zip(samples['parent_name'], samples['reference_file']):
        parents[k] = v
    
    # if one of the parents, return parent sequences
    if wildcards.sample in parents.keys():
        return parents[wildcards.sample]
    
    # otherwise, just return reference
    return get_column_by_sample(wildcards, samples, 'reference_file')

# map to one of the parental references.  The choice of reference is arbitrary
rule align:
    input:
        reads = get_reads,
        reference = get_reference
    output:
        aligned = "out/aligned/{sample}.bam",
        idx = "out/aligned/{sample}.bam.bai",
    conda: "../deps/align/env.yml"
    container: "docker://szsctt/lr_align"
    threads: 8
    shell:
        """
        minimap2 -t {threads} -ax map-hifi {input.reference} {input.reads} -B 1.5 --end-bonus 5 --MD |\
            samtools sort -o {output.aligned} -
        
        samtools index {output.aligned}
        """