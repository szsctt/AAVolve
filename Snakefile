from scripts.get_samples import get_samples

samples = get_samples(config)

wildcard_constraints:
    sample_name = '|'.join(samples.sample_name),



consensus = list()
for name, seq_tech in zip(samples.sample_name, samples.seq_tech):
    if seq_tech == 'np-cc':
        consensus.append(f"out/c3poa_filt/{name}.fasta.gz")
        consensus.append(f"out/c3poa/{name}/repeat_counts.tsv")

rule all:
    input: consensus

include: 'rules/consensus.smk'