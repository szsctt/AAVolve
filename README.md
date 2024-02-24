# AAVolve
AAVolve: A tool for exploring directed evolution of AAV capsids using long-read sequencing data

AAVolve is a snakemake pipeline for analysing long-read, shuffled AAV capsid data.  It supports nanopore (R2C2) and PacBio data.


## Quickstart

To run with a local installation of snakemake and singularity, use:

```
snakemake --cores 1 --config read_file=<path to fastq> parent_file=<path to fasta> seq_tech=<np, np-cc, or pb>
```

For Nanopore R2C2 data (np-cc), please also specify a path to a fasta file containing the splint:

```
snakemake --cores 1 --config read_file=<path to fastq> parent_file=<path to fasta> seq_tech=<np, np-cc, or pb> splint_file=<path_to_fasta>
```

## Requirements

For more than one sample, specify paramters and inputs in a comma-sepearated file with one row per sample (fastq file), and the following columns:

sample_name,parent_name,reference_name,seq_tech,min_reps,read_file,parent_file,reference_file,splint_file,non_parental_freq

- `sample_name`: A name for the sample
- `parent_name` (optional): A name for the parents of the sample
- `reference_name` (optional): A name for the parent to be used as the reference
- `seq_tech`: The sequencing technology used: either `np` for Nanopore, `np-cc` for Nanopore RCA/R2C2, or `pb` for PacBio
- `min_reps` (optional): For RCA data, the minimum number of repeats for a read to be included (default: 3)
- `read_file`: The path to the fastq (or fasta) file containing the data for the sample
- `parent_file`: The path to the fasta file containing the parental sequences used for the shuffling experiment
- `reference_file` (optional): The path to the fasta file containing one of the parental sequences to be used as a reference for alingment (default: the first parental sequence)
- `splint_file` (required for RCA data): The path to the fasta file containing the splint used for circularization during RCA
- `non_parental_freq` (optinal): Include non-parental variants seen in more than this fraction of the reads (default: don't include non-parental variants)


