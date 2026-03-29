import os
import re

from aavolve.get_samples import get_samples
from aavolve.snakemake_helpers import (
    build_group_report_targets,
    build_input_validation_targets,
    build_npcc_consensus_maps,
    build_non_parental_variant_group_maps,
)
from aavolve.utils import normalize_names


samples = get_samples(config)

# When running in single-sample mode (no samples CSV was provided), write one so
# the group_manifest rule (which requires a CSV path) can still function.
if 'samples' not in config:
    _samples_csv = os.path.join("out", "generated_samples.csv")
    os.makedirs(os.path.dirname(_samples_csv), exist_ok=True)
    samples.to_csv(_samples_csv, index=False)
    config['samples'] = _samples_csv

sample_names = normalize_names(samples.sample_name)
parent_names = normalize_names(samples.parent_name)
all_names = sorted(set(sample_names) | set(parent_names))
if not all_names:
    raise ValueError("No sample or parent names available to constrain the 'sample' wildcard")

sample_name_pattern = "|".join(re.escape(name) for name in all_names)

wildcard_constraints:
    sample=sample_name_pattern,
    anchor_suffix="\.fastq(\.gz)?|\.fasta(\.gz)?",
    # Derived key used to de-duplicate expensive np-cc (R2C2) consensus across samples.
    cc_id="npcc_[0-9a-f]{12}",

# target files for RCA consensus
consensus = list()
for name, seq_tech in zip(samples.sample_name, samples.seq_tech):
    if seq_tech == 'np-cc':
        consensus.append(f"out/c3poa_filt/{name}.fasta.gz")
        consensus.append(f"out/c3poa/{name}/repeat_counts.tsv")

npcc_sample_to_cc_id, npcc_cc_id_to_reads, npcc_cc_id_to_splint = build_npcc_consensus_maps(samples)


# Validate inputs once per unique (parent_file, reference_file) combination.
input_validation_map, input_validation_samples, input_validation_targets = build_input_validation_targets(samples)
_, _, group_report_targets = build_group_report_targets(samples)

sample_to_non_parental_group_id, non_parental_variant_groups = build_non_parental_variant_group_maps(
    samples, input_validation_samples
)
non_parental_variant_group_ids = sorted(
    group_id for group_id, group_samples in non_parental_variant_groups.items() if len(group_samples) > 1
)

trim_qc_targets = []
if "trim" in samples.columns:
    trim_qc_targets = [
        f"out/qc/trimming/{sample}.mafft.fasta"
        for sample, trim in zip(samples.sample_name, samples.trim)
        if isinstance(trim, bool) and trim
    ]

pretrim_qc_targets = [f"out/qc/trimming/{sample}.pretrim.mafft.fasta" for sample in samples.sample_name]

rule all:
    input: 
        consensus,
        input_validation_targets,
        group_report_targets,
        pretrim_qc_targets,
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
