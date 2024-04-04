import numpy as np
from snakemake.io import expand

# general helpers
def get_column_by_sample(wildcards, samples, column_name):

    assert len(samples.sample_name) == len(samples.sample_name.unique()), "Sample names are not unique"
    return {k:v for k, v in zip(samples.sample_name, samples[column_name])}[wildcards.sample]

def is_fastq(file):
    return any((file.endswith('.fastq'), file.endswith('.fastq.gz'), file.endswith('.fq'), file.endswith('.fq.gz')))


#### align ####

def get_reads(wildcards, samples):
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
    if tech == 'np-cc':
        
        # check if we want to filter for repeats
        repeats = get_column_by_sample(wildcards, samples, 'min_reps')

        if np.isnan(repeats):
            return f"out/c3poa/{wildcards.sample}/split/R2C2_Consensus.fasta.gz"

        else:
            return f"out/c3poa_filt/{wildcards.sample}.fasta.gz"

    # otherwise, just return reads
    return get_column_by_sample(wildcards, samples, 'read_file')

def get_reference(wildcards, samples):
    """
    Get appropriate reference for wildcards.sample
    Either parental references,
    or just the reference otherwise
    """
    if wildcards.sample not in set(samples.sample_name) | set(samples.parent_name):
        raise ValueError(f"Sample {wildcards.sample} not found")

    # make a dictionary of parents
    parents = {}
    for k, v in zip(samples['parent_name'], samples['reference_file']):
        parents[k] = v
    
    # if one of the parents, return parent sequences
    if wildcards.sample in parents.keys():
        return parents[wildcards.sample]
    
    # otherwise, just return reference
    return get_column_by_sample(wildcards, samples, 'reference_file')

#### variants ####

# rules num_parents
def get_parents(wildcards, samples):
    
    if wildcards.sample not in set(samples.parent_name):
        raise ValueError(f"Parent {wildcards.sample} not found")

    parents = {k:v for k,v in zip(samples.parent_name, samples.parent_file)}[wildcards.sample]
    
    # check there is a value for the parent file
    if parents is None:
        raise ValueError(f"Sample {wildcards.sample} does not have a parent")
    try:
        if np.isnan(parents):
            raise ValueError(f"Sample {wildcards.sample} does not have a parent")
    except TypeError:
        pass

    return parents


def fill_parents(wildcards, samples, filename):
    parent_name = get_column_by_sample(wildcards, samples, "parent_name")
    return expand(filename, sample=parent_name)

#### tranform_variants ###

# rule dmat
def get_dmat_input(wildcards, nt_seq_file, aa_seq_file):
    if wildcards.seq_type == "nt-seq":
        return nt_seq_file
    elif wildcards.seq_type == "aa-seq":
        return aa_seq_file
    else:
        raise ValueError("seq_type must be nt-seq or aa-seq")

# rule count_reads
def get_reads_for_counting(wildcards, samples, consensus, con_filt):
    """
    Get the appropriate files for counting original reads
    If not R2C2 data, just return result of get_reads(wildcards)

    If R2C2 data, return original reads and consensus reads
    If we filtered, also return filtered reads
    """

    # get input reads
    reads = [{k:v for k, v in zip(samples.sample_name, samples.read_file)}[wildcards.sample]]
    
    # get sequencing technology
    tech = {k:v for k, v in zip(samples.sample_name, samples.seq_tech)}[wildcards.sample]
    # if not R2C2, just return reads
    if tech != 'np-cc':
        return reads

    # for R2C2 add consensus reads
    reads.append(consensus)
    # and filtered consensus reads if they exist
    filt = {k:v for k, v in zip(samples.sample_name, samples.min_reps)}[wildcards.sample]
    if filt is not None and np.isnan(filt) == False:
        reads.append(con_filt)
    
    return reads

# rule count_reads
def format_input_reads(input):
    """
    Format input reads for counting
    """
    fastqs = []
    fastas = []
    
    # check if input was fastq (most reads) or fasta (sanger sequencing)
    for i in input:
        if is_fastq(i):
            fastqs.append(i)
        else:
            fastas.append(i)
    
    # format output
    args = ''
    if len(fastas) > 0:
        args = args + f"--fasta-files {' '.join(fastas)} "
    if len(fastqs) > 0:
        args = args + f"--fastq-files {' '.join(fastqs)}"

    return args