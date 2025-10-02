# TideHunter vs C3POa Comparison

This directory contains scripts and outputs for comparing **TideHunter** and **C3POa** performance on R2C2 (Rolling Circle Amplification to Concatemeric Consensus) reads from AAV2.

## Purpose

Compare the two consensus-calling tools for processing circular consensus sequencing data:
- **TideHunter**: A tool for identifying tandem repeats in long reads
- **C3POa**: Concatemeric Consensus Caller with Partial Order Alignments, specifically designed for R2C2 data

## Test Data

- **Source**: R2C2 reads from AAV2 (Adeno-Associated Virus serotype 2)
- **Technology**: Nanopore long-read sequencing with Rolling Circle Amplification
- **Reference**: AAV2 capsid sequence (2,208 bp)

## Directory Structure

```
tidehunter_test/
├── data/                      # Input data and references
├── output/                    # Generated outputs
│   ├── tidehunter/           # TideHunter consensus sequences
│   ├── c3poa/                # C3POa consensus sequences
│   ├── alignments/           # BAM files (aligned to AAV2 reference)
│   └── plots/                # Comparison plots and statistics
└── scripts/                  # Analysis scripts
    ├── run.sh                # Master pipeline script
    ├── compute_accuracy_by_repeats.py
    └── plot_repeats.py
```

## Analysis Pipeline

The `run.sh` script orchestrates the full analysis:

1. **Environment Setup**
   - Creates/activates `tidehunter_env` micromamba environment
   - Installs required dependencies (matplotlib, pysam, etc.)

2. **Consensus Calling**
   - **TideHunter** (two modes):
     - With primers: `-5 PRIMER_SEQ -3 PRIMER_SEQ`
     - Without primers: no arguments (identifies repeats empirically)
   - **C3POa**: Uses splint sequence for R2C2 processing

3. **Alignment**
   - Aligns consensus sequences to AAV2 reference using minimap2
   - Generates sorted, indexed BAM files

4. **Comparative Analysis**
   - Repeat count distributions
   - Accuracy (percent identity) by number of repeats
   - Aligned length by number of repeats
   - Full-reference coverage by number of repeats

## Key Scripts

### `compute_accuracy_by_repeats.py`

Analyzes BAM alignments and computes per-repeat statistics:

**Usage:**
```bash
python3 compute_accuracy_by_repeats.py \
  --tide-bam tidehunter_empiricalseqs.bam \
  --tide-bam-noargs tidehunter_noprimers.bam \
  --c3poa-bam c3poa.bam \
  --out-tsv accuracy_by_repeats.tsv \
  --plot accuracy_by_repeats.png \
  --length-plot length_by_repeats.png \
  --cover-plot cover_by_repeats.png \
  --log-x
```

**Features:**
- Extracts repeat counts from read headers
- Computes percent identity using CIGAR strings and NM tags
- Defines full-reference coverage using alignment span (reference_start ≤ 0 and reference_end ≥ ref_len)
- Adds +1 to C3POa repeat counts to match TideHunter convention (C3POa uses 0 for single-pass, TideHunter uses 1)
- Supports log-scale x-axis with automatic shifting for zero/negative values
- Adds vertical reference line at 3 repeats on all plots

**Outputs:**
- TSV file with columns: `tool`, `repeats`, `n_reads`, `mean_id`, `median_id`, `sd_id`, `median_aligned_len`, `frac_full_cover`
- PNG plots:
  - Median percent identity vs repeats
  - Median aligned length vs repeats
  - Number of reads covering full reference vs repeats
- Summary statistics: Total reads with ≥3 repeats covering full reference

### `plot_repeats.py`

Generates frequency polygon comparing repeat count distributions:

**Usage:**
```bash
python3 plot_repeats.py \
  --tidehunter tidehunter_output.fa \
  --tidehunter-noargs tidehunter_noprimers_output.fa \
  --c3poa R2C2_Consensus.fasta \
  --out repeats_compare.png \
  --bins 80 \
  --log-x
```

**Features:**
- Parses repeat counts from FASTA headers
- Adds +1 to C3POa counts for consistency
- Creates frequency polygon (line plot through histogram bin centers)
- Optional log-scale x-axis
- Vertical reference line at 3 repeats

## Key Findings Format

The analysis produces:

1. **Repeat Distribution Plot** (`repeats_compare.png`)
   - Shows how many reads have each repeat count for each tool

2. **Accuracy Analysis** (`accuracy_by_repeats.png`)
   - Median percent identity for each repeat count
   - Higher repeats generally mean higher accuracy

3. **Length Analysis** (`length_by_repeats.png`)
   - Median aligned length by repeat count
   - Shows how alignment length correlates with repeat count

4. **Coverage Analysis** (`cover_by_repeats.png`)
   - Number of reads covering full reference (start ≤ 0, end ≥ 2208 bp)
   - Critical for downstream analysis requiring complete sequences

5. **Summary Statistics**
   ```
   === Summary: Reads with >=3 repeats covering full reference ===
   tidehunter:         89,415 reads
   tidehunter_noargs:  9,196 reads
   c3poa:              136,573 reads
   ```

## Important Implementation Details

### Repeat Count Convention
- **TideHunter**: Uses 1 for single-pass reads, 2 for one repeat, etc.
- **C3POa/R2C2**: Uses 0 for single-pass reads, 1 for one repeat, etc.
- **Solution**: Add +1 to all C3POa repeat counts for direct comparison

### Coverage Definition
Full-reference coverage is defined as:
- `alignment.reference_start ≤ 0` AND `alignment.reference_end ≥ reference_length`
- Fallback to aligned base count (≥ ref_len) if coordinates unavailable
- More stringent than just counting aligned bases

### Percent Identity Calculation
Using pysam's `get_cigar_stats()`:
```python
aligned_bases = M_count + match_count + mismatch_count  # indices 0, 7, 8
matches = aligned_bases - NM_tag  # NM is edit distance
percent_identity = (matches / aligned_bases) * 100
```

### Log-Scale Handling
When `--log-x` is used:
- If `min(repeats) ≤ 0`, shift all values by `1 - min(repeats)`
- Apply `plt.xscale('log')`
- Adjust vertical reference line position accordingly
- Label x-axis with shift amount if applied

## Running the Full Pipeline

```bash
# From the AAVolve repository root
cd tidehunter_test
bash scripts/run.sh
```

The script will:
1. Set up the environment (if needed)
2. Run TideHunter (both modes)
3. Run C3POa
4. Align all outputs to AAV2 reference
5. Generate comparative plots and statistics

## Dependencies

- **Python 3.10+** with:
  - matplotlib
  - numpy
  - pysam
  - biopython (for aavolve.utils)
  
- **Bioinformatics tools**:
  - minimap2
  - samtools
  - TideHunter
  - C3POa

- **Environment manager**: micromamba or conda

## References

- **TideHunter**: Gao et al. (2019) "TideHunter: efficient and sensitive tandem repeat detection from noisy long reads"
- **C3POa**: Volden et al. (2018) "Improving nanopore read accuracy with the R2C2 method enables the sequencing of highly multiplexed full-length single-cell cDNA"
- **R2C2**: Rolling Circle Amplification to Concatemeric Consensus sequencing method

## Notes

- All plots include a dotted vertical line at 3 repeats as a visual reference
- The pipeline uses alignment span (not just aligned bases) to define full-reference coverage
- TideHunter modes:
  - With primers: More accurate when primer sequences are known
  - Without primers: Useful when primer sequences are unknown or variable
- C3POa requires a splint sequence for proper R2C2 processing
