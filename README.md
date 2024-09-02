# AAVolve
AAVolve: A tool for exploring directed evolution of AAV capsids using long-read sequencing data

AAVolve is a snakemake pipeline for analysing long-read, shuffled AAV capsid data.  It supports Nanopore (including R2C2), PacBio (including HiFi) and Sanger data.

## Quickstart

To run with a local installation of `snakemake` and `apptainer`, clone this repository and use:

```
snakemake --use-apptainer --cores 1 --config read_file=<path to fastq> parent_file=<path to fasta> seq_tech=<np, np-cc, pb or sg>
```

For Nanopore R2C2 data (np-cc), please also specify a path to a fasta file containing the splint:

```
snakemake --use-apptainer --cores 1 --config read_file=<path to fastq> parent_file=<path to fasta> seq_tech=<np, np-cc, or pb> splint_file=<path_to_fasta>
```

Look for outputs in the `out` folder.

## Setup and usage

AAVolve is a [`snakemake`](https://snakemake.readthedocs.io/en/stable/) pipeline.  The easist way to use AAVolve is to tell `snakemake` to use [`apptainer`](https://apptainer.org/) (previously known as singularity) to automatically supply dependencies. For this, please first install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (we recommend using [conda/mamba](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba)] and [`apptainer`](https://apptainer.org/docs/user/latest/quick_start.html#installation).  Note that Apptainer is only available for Linux.

With `snakemake` and `apptainer` installed, clone this repository, and change into the cloned directory: 

```
git clone https://github.com/szsctt/AAVolve.git
cd AAVolve
```

AAVolve makes use of `snakemake's` `--config` or `--configfile` argument to specify the input files and parameters for analysis. At a minimum, specify a file containing reads (`read_file`; fgastq format), a file containg parental sequences (`parent_file`; fasta format), and the sequencing technology used (`seq_tech`; `np` for Nanopore, `np-cc` for Nanopore R2C2, `pb` for PacBio and `sg` for Sanger). For Nanopore R2C2 data, a file containing the splint sequence (`splint_file`; fasta format) must also be provided. Other optional config parameters can be specified as outlined below ('Config options'). To analyse one file, use the `--config` parameter to specify inputs on the command line; to analyse multiple files it may be more convenient to specify inputs in a comma-separated file as outlined in the section below ('Config options').

`snakemake` also requires the user to specify the number of cores to be used with the `--cores` argument.

To run AAVolve using one core, for one read file and using apptainer to supply dependencies:

```
snakemake --use-apptainer --cores 1 --config read_file=<path to fastq> parent_file=<path to fasta> seq_tech=<np, np-cc, or pb> splint_file=<path_to_fasta>
```

Where the `splint_file` argument can be omitted for non Nanopore R2C2 data.

See the [`snakemake` documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more command line options, including cluster execution and profiles.

### Config options

For more than one sample, specify parameters and inputs in a comma-sepearated file with one row per sample (fastq file), and the following columns:

sample_name,parent_name,reference_name,seq_tech,min_reps,read_file,parent_file,reference_file,splint_file,non_parental_freq

- `sample_name`: A name for the sample
- `parent_name` (optional): A name for the parents of the sample
- `reference_name` (optional): A name for the parent to be used as the reference
- `seq_tech`: The sequencing technology used: either `np` for Nanopore, `np-cc` for Nanopore RCA/R2C2, `pb` or `pb-hifi` for PacBio, or `sg` for sanger sequencing (fasta format). 
- `min_reps` (optional): For RCA data, the minimum number of repeats for a read to be included (default: 3)
- `read_file`: The path to the fastq (or fasta) file containing the data for the sample. For PacBio HiFi data, provide the consensus (ccs) reads.
- `parent_file`: The path to the fasta file containing the parental sequences used for the shuffling experiment
- `reference_file` (optional): The path to the fasta file containing one of the parental sequences to be used as a reference for alingment (default: the first parental sequence)
- `splint_file` (required for RCA data): The path to the fasta file containing the splint used for circularization during RCA
- `non_parental_freq` (optional): Identify non-parental variants seen in more than this fraction of the reads (default: 0.2)
- `include_non_parental`(optional): Include non-parental variants if they are seen in more than `non_parental_freq` reads (default: False)
- `group_vars` (optional): During assignment of parents, group variants and assign parent with lowest hamming distance to read.  If False, variants are considered individually (default: True)
- `group_vars_dist` (optional): When grouping variants, combine non-adjacent variants that are at most this far apart (defult: 4)
- `max_group_distance` (optional): When grouping variants, assign a parent if there are at most this fraction of the variants that differ between a parent and the read. For example, in a group of 6 variants, if this parameter is set to 0.2, a parent that has one variant that differs from the read will still be assigned to the group.  Note that if `group_vars` is `False`, variants must match a parent otherwise they will be set to `NA`, regardless of the value for this parameter (default: 0.2).

To provide this file to `snakemake`:

```
snakemake --use-apptainer --cores 1 --config samples=<path to csv>
```

## Outputs

AAVolve produces several outputs, which will appear in the `out` folder in  may be useful depending on the library to be analysed:

- Summary of results
- 

