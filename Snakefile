from aavolve.get_samples import get_samples

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
        expand("out/variants/reads/{sample}.tsv.gz", sample=samples.sample_name),
        expand("out/variants/reads/{sample}_read-count.txt", sample=samples.sample_name),
        expand("out/variants/combined/{sample}.tsv.gz", sample=samples.sample_name),
        expand("out/variants/pivot/{sample}_parents.tsv.gz", sample=samples.sample_name),
        expand("out/parents/assigned/{sample}_assigned-parents.tsv.gz", sample=samples.sample_name),
        expand("out/parents/freqs/{sample}_assigned-parents_freq.tsv.gz", sample=samples.sample_name),
        expand("out/parents/breaks/{sample}.tsv.gz", sample=samples.sample_name),
        expand("out/corrected/counts/{sample}_{seqtype}-seq-counts.tsv.gz", sample=samples.sample_name, seqtype = ("aa", "nt")),
        expand("out/corrected/dmat/{sample}_{subset}_{seqtype}-seq.tsv.gz", sample=samples.sample_name, subset = ("random", "first"), seqtype = ("aa", "nt")),
        expand("out/qc/{sample}_read-counts.tsv", sample=samples.sample_name),
        expand("out/qc/{sample}_report.html", sample=samples.sample_name),

include: 'rules/consensus.smk'
include: 'rules/align.smk'
include: 'rules/variants.smk'
include: 'rules/transform_variants.smk'
