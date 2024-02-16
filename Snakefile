from scripts.get_samples import get_samples

samples = get_samples(config)

wildcard_constraints:
    sample_name = '|'.join(samples.sample_name) + "|" + "|".join(samples.parent_name)

# target files for RCA consensus
consensus = list()
for name, seq_tech in zip(samples.sample_name, samples.seq_tech):
    if seq_tech == 'np-cc':
        consensus.append(f"out/c3poa_filt/{name}.fasta.gz")
        consensus.append(f"out/c3poa/{name}/repeat_counts.tsv")

rule all:
    input: 
        consensus,
        expand("out/aligned/{sample}.bam", sample=samples.sample_name),
        expand("out/aligned/{sample}.bam", sample=samples.parent_name),
        expand("out/variants_parents/{sample}.tsv.gz", sample=samples.parent_name)

include: 'rules/consensus.smk'
include: 'rules/align.smk'
include: 'rules/variants.smk'
