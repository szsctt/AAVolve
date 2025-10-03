# Adapter Trimming Comparison Test

This directory contains a benchmark comparison of different adapter trimming tools for long-read sequencing data from AAV (Adeno-Associated Virus) libraries sequenced using Nanopore R2C2 (Rolling Circle to Concatemeric Consensus) technology.

## Overview

The test evaluates the performance and effectiveness of multiple adapter trimming tools on consensus reads generated from R2C2 sequencing data. The goal is to identify the most suitable tool for removing splint/adapter sequences from AAV capsid library sequences.

## Dataset

- **Source**: SRA accession SRR29543177
- **Sample**: AAV2 library with R2C2 sequencing
- **Technology**: Nanopore R2C2/RCA (Rolling Circle Amplification)
- **Preprocessing**: C3POa consensus calling to generate high-quality consensus reads
- **Total reads**: 275,458 consensus reads after C3POa processing

### Adapter Sequences

The splint sequences used for R2C2 sequencing that need to be trimmed:
- **5' adapter**: `TGATTTAAATCAGGT` (15 bp)
- **3' adapter**: `TTGCTTGTTAAATCA` (15 bp)

These sequences are stored in `data/adapters.fa`.

## Tools Tested

### 1. Trimmomatic
- **Mode**: Single-end (SE) with ILLUMINACLIP
- **Parameters**: `ILLUMINACLIP:adapters.fa:2:30:10`
- **Characteristics**: 
  - Originally designed for Illumina data
  - Reports reads dropped vs. surviving, but not trimming statistics
  - Does not report how many reads were actually trimmed

### 2. Cutadapt (Internal mode)
- **Mode**: Internal adapter search (default)
- **Parameters**: `-a TTGCTTGTTAAATCA -g TGATTTAAATCAGGT --revcomp --trimmed-only`
- **Characteristics**:
  - Searches for adapters anywhere within reads
  - Handles reverse complements automatically
  - Reports detailed trimming statistics

### 3. Cutadapt (Non-internal mode)
- **Mode**: Non-internal adapter search (anchored to ends)
- **Parameters**: `-a TTGCTTGTTAAATCAX -g XTGATTTAAATCAGGT --revcomp --trimmed-only`
- **Characteristics**:
  - Only trims adapters at read ends (X prefix indicates non-internal)
  - More conservative approach
  - Faster than internal mode

### 4. Cutadapt (Linked adapters)
- **Mode**: Linked adapter pairs
- **Parameters**: `-a TGATTTAAATCAGGT...TTGCTTGTTAAATCA --revcomp --trimmed-only`
- **Characteristics**:
  - Requires both adapters to be present in order
  - Useful for RCA data where splints should flank the insert
  - Most stringent filtering

### 5. BBDuk (Left trim)
- **Mode**: Left k-mer trimming
- **Parameters**: `ktrim=l k=13 mink=5 hdist=1 rcomp=t`
- **Characteristics**:
  - K-mer based matching (k=13, minimum k=5)
  - Allows 1 mismatch (hdist=1)
  - Trims from left side of reads

### 6. BBDuk (Right trim)
- **Mode**: Right k-mer trimming (applied after left trim)
- **Parameters**: `ktrim=r k=13 mink=5 hdist=1 rcomp=t`
- **Characteristics**:
  - Applied sequentially after left trim
  - Same k-mer matching parameters
  - Trims from right side of reads

## Directory Structure

```
trimming_test/
├── README.md                    # This file
├── data/                        # Input data and references
│   ├── SRR29543177.fastq       # Raw R2C2 reads (downloaded from SRA)
│   ├── aav2.fa                 # AAV2 reference sequence
│   ├── adapters.fa             # Adapter sequences for trimming
│   └── splint.fa               # Splint sequence for C3POa
├── scripts/
│   └── run.sh                  # Main pipeline script
└── output/                      # Results for each tool
    ├── c3poa/                  # C3POa consensus reads
    ├── trimmomatic/            # Trimmomatic results
    ├── cutadapt/               # Cutadapt internal mode results
    ├── cutadapt_noninternal/   # Cutadapt non-internal mode results
    ├── cutadapt_linked/        # Cutadapt linked adapters results
    ├── bbmap/                  # BBDuk results (both left and right trim)
    └── msa/                    # Multiple sequence alignments (first 100 reads)
```

## Running the Pipeline

### Prerequisites

- **Singularity/Apptainer**: Required for C3POa
- **Micromamba/Conda**: For managing software dependencies
- Run from the **project root directory** (AAVolve/)

### Quick Start

```bash
# From the AAVolve project root:
bash trimming_test/scripts/run.sh
```

The script will:
1. Create a conda environment with required tools
2. Download the test dataset from SRA (if not present)
3. Run C3POa to generate consensus reads (currently commented out)
4. Convert FASTA to FASTQ (Trimmomatic requirement)
5. Run all trimming tools with timing measurements
6. Generate multiple sequence alignments for visual comparison
7. Display a summary of results

### Environment Setup

The script automatically creates a conda environment with:
- `sra-tools` (for downloading data)
- `mafft` (for multiple sequence alignment)
- `cutadapt` (adapter trimming)
- `trimmomatic` (adapter trimming)
- `bbmap` (BBDuk for adapter trimming)

## Output Files

Each tool's output directory contains:
- **`*_trimmed.fastq`** or **`*_cutadapt.fastq`**: Trimmed reads
- **`*_time.txt`**: Runtime statistics and tool-specific output from `/usr/bin/time -v`

The `msa/` directory contains:
- **`{tool}_first100.fa`**: First 100 reads for each tool (temporary)
- **`{tool}_with_aav2.fa`**: Reads concatenated with AAV2 reference
- **`{tool}_aligned.fa`**: MAFFT multiple sequence alignment
- **`{tool}_mafft.log`**: MAFFT alignment log

## Results Summary

The script generates a summary showing for each tool:
- **Runtime**: Wall clock time for processing
- **Total reads**: Number of input reads
- **Surviving reads**: Reads retained after processing
- **Reads trimmed/dropped**: Number of reads modified or removed

### Example Output

```
========================================
TRIMMING SUMMARY
========================================

Trimmomatic:
  Runtime: 0:05.27
  Total reads: 275458
  Surviving reads: 275458
  Dropped reads: 0

Cutadapt (internal):
  Runtime: 0:47.42
  Total reads: 275458
  Surviving reads: 274286
  Reads with adapters trimmed: 274286

Cutadapt (non-internal):
  Runtime: 0:02.37
  Total reads: 275458
  Surviving reads: 456
  Reads with adapters trimmed: 456

Cutadapt (linked):
  Runtime: 0:47.81
  Total reads: 275458
  Surviving reads: 274286
  Reads with adapters trimmed: 274286

BBDuk (left trim):
  Runtime: 0:03.70
  Total reads: 275458
  Surviving reads: 275458
  Reads trimmed: 269968

BBDuk (right trim):
  Runtime: 0:02.73
  Total reads: 275458
  Surviving reads: 275458
  Reads trimmed: 10406

========================================
```

## Key Findings

1. **Speed**: BBDuk is fastest (2.7-3.7s), followed by Cutadapt non-internal (2.4s) and Trimmomatic (5.3s). Cutadapt internal and linked modes are slowest (~47s).

2. **Adapter Detection and Filtering** (with `--trimmed-only` flag):
   - **Cutadapt internal**: Found adapters in 274,286/275,458 reads (99.6%), filtered out 1,172 reads without adapters
   - **Cutadapt non-internal**: Only found adapters at read ends in 456 reads (0.17%), filtered out 274,002 reads
   - **Cutadapt linked**: Same detection as internal mode (274,286 reads, 99.6%)
   - **BBDuk left trim**: Trimmed 269,968 reads (98.0%), kept all reads
   - **BBDuk right trim**: Trimmed only 10,406 reads (3.8%) after left trim, kept all reads
   - **Trimmomatic**: Kept all 275,458 reads (no filtering)

3. **Impact of `--trimmed-only` Flag**: 
   - This flag causes Cutadapt to discard reads without detected adapters
   - Explains why Cutadapt outputs have fewer surviving reads than input reads
   - Non-internal mode is most aggressive, keeping only 0.17% of reads (those with adapters strictly at ends)
   - Internal and linked modes keep 99.6% of reads (those with adapters found anywhere)

4. **Sequential BBDuk Trimming**:
   - Left trim removes adapters from 98% of reads
   - Right trim (applied after left trim) only finds adapters in 3.8% of reads
   - Suggests most adapters are on the left side, or that left trim removes most adapter sequences
   - Both steps retain all reads (no filtering)

5. **Tool Reporting Differences**: 
   - **Cutadapt**: Reports detailed adapter statistics and base pair changes
   - **BBDuk**: Reports k-mer trimming statistics and base pair changes
   - **Trimmomatic**: Only reports survival/dropout counts, not how many reads were actually trimmed

## Notes

- FASTQ conversion is required because Trimmomatic only accepts FASTQ format (dummy Q40 quality scores added)
- Multiple sequence alignments use only the first 100 reads for performance and visualization purposes
- The `--trimmed-only` flag in Cutadapt means only reads with adapters found are written to output

## References

- **C3POa**: [https://github.com/rvolden/C3POa](https://github.com/rvolden/C3POa)
- **Cutadapt**: [https://cutadapt.readthedocs.io/](https://cutadapt.readthedocs.io/)
- **Trimmomatic**: [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- **BBDuk (BBMap suite)**: [https://jgi.doe.gov/data-and-tools/software-tools/bbtools/](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
- **Dataset**: SRA SRR29543177

