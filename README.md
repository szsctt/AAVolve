# AAVolve
AAVolve: A tool for exploring directed evolution of AAV capsids using long-read sequencing data

AAVolve is a snakemake pipeline for analysing long-read, shuffled AAV capsid data.  It supports nanopore (R2C2) and PacBio data.


## Quickstart

To run with a local installation of snakemake and singularity, use:

```
snakemake --cores 1 --config read_file=<path to fastq> parent_file=<path to fasta> seq_tech=<np, np-cc, or pb>
```

## Requirements




