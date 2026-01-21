import re

from aavolve.get_samples import get_samples
from aavolve.snakemake_helpers import build_group_report_targets, build_input_validation_targets
from aavolve.utils import normalize_names


samples = get_samples(config)

sample_names = normalize_names(samples.sample_name)
parent_names = normalize_names(samples.parent_name)
all_names = sorted(set(sample_names) | set(parent_names))
if not all_names:
    raise ValueError("No sample or parent names available to constrain the 'sample' wildcard")

sample_name_pattern = "|".join(re.escape(name) for name in all_names)

wildcard_constraints:
    sample=sample_name_pattern,
    anchor_suffix="\.fastq(\.gz)?|\.fasta(\.gz)?",

# target files for RCA consensus
consensus = list()
for name, seq_tech in zip(samples.sample_name, samples.seq_tech):
    if seq_tech == 'np-cc':
        consensus.append(f"out/c3poa_filt/{name}.fasta.gz")
        consensus.append(f"out/c3poa/{name}/repeat_counts.tsv")


# Validate inputs once per unique (parent_file, reference_file) combination.
input_validation_map, input_validation_samples, input_validation_targets = build_input_validation_targets(samples)
_, _, group_report_targets = build_group_report_targets(samples)

trim_qc_targets = []
if "trim" in samples.columns:
    trim_qc_targets = [
        f"out/qc/trimming/{sample}.mafft.fasta"
        for sample, trim in zip(samples.sample_name, samples.trim)
        if isinstance(trim, bool) and trim
    ]

rule all:
    input: 
        consensus,
        input_validation_targets,
        group_report_targets,
        trim_qc_targets,
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
        expand("out/reports/sample_reports/{sample}_report.html", sample=samples.sample_name),

include: 'rules/consensus.smk'
include: 'rules/check_inputs.smk'
include: 'rules/align.smk'
include: 'rules/variants.smk'
include: 'rules/transform_variants.smk'
include: 'rules/group_report.smk'
