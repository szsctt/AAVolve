# AAVolve Copilot Instructions

AAVolve is a Snakemake pipeline for analyzing long-read sequencing data from directed evolution experiments with shuffled AAV (Adeno-Associated Virus) capsid libraries.


CORRECTNESS IS VERY IMPORTANT. Do not introduce fallbacks, it's better to fail than to produce incorrect results.

## Project Architecture

### Core Components
- **Snakemake Pipeline**: Main workflow orchestrator in `Snakefile` and `rules/` directory
- **Python Package**: `aavolve/` contains modular processing scripts for bioinformatics operations
- **Containerization**: Uses Docker/Apptainer containers (docker://szsctt/lr_pybio:py310) for reproducible dependencies
- **Config-driven**: CSV configuration files define samples, sequencing technologies, and parameters

### Key Data Flow
1. **Input**: FASTQ/FASTA files with parent sequences (AAV capsid variants)
2. **Consensus** (RCA only): C3POa consensus calling for np-cc (Nanopore RCA/R2C2) data
3. **Alignment**: minimap2 alignment against reference sequences 
4. **Variant Calling**: Extract variants from SAM files, compare to parental sequences
5. **Parent Assignment**: Assign reads to most likely parental origin based on variant patterns
6. **Analysis**: Count unique sequences, generate distance matrices, frequency analysis

### Sequencing Technology Support
- `np`: Nanopore
- `np-cc`: Nanopore RCA/R2C2 (requires splint sequence)
- `pb`/`pb-hifi`: PacBio HiFi
- `sg`: Sanger sequencing

## Development Patterns

### Snakemake Rules Structure
```python
# Rules are organized in separate files under rules/
include: 'rules/consensus.smk'    # C3POa consensus for RCA data
include: 'rules/align.smk'        # Minimap2 alignment
include: 'rules/variants.smk'     # Variant extraction and calling
include: 'rules/transform_variants.smk'  # Downstream analysis
```

### Python Module Pattern
Each `aavolve/*.py` script follows the pattern:
- Command-line interface using argparse
- Main function for programmatic use
- Can be called as `python3 -m aavolve.script_name`
- Example: `python3 -m aavolve.extract_features_from_sam -i input.bam -r ref.fa -o output.tsv`

### Configuration System
- Use `aavolve.get_samples.get_samples(config)` to parse configuration
- Supports both command-line `--config` and CSV `--configfile` modes
- Required fields: `sample_name`, `parent_file`, `read_file`, `seq_tech`
- Wildcard constraints generated from sample names to ensure proper rule matching

### Testing Conventions
- pytest framework in `tests/pytest/`
- One test file per Python module (e.g., `test_extract_features_from_sam.py`)
- Use `conftest.py` fixtures for sample data and configurations
- Test data in `tests/data/` directory
- Run tests: `pytest tests/pytest/`

## Common Tasks

### Adding New Analysis Module
1. Create `aavolve/new_module.py` with main() function and argparse CLI
2. Add corresponding `test_new_module.py` in `tests/pytest/`
3. Add rule to appropriate `rules/*.smk` file if part of pipeline
4. Update wildcard constraints if new sample types introduced

### Running Pipeline
```bash
# Single sample via command line
snakemake --use-apptainer --cores 1 --config read_file=data.fastq parent_file=parents.fa seq_tech=np sample_name=sample1

# Multiple samples via CSV config
snakemake --use-apptainer --cores 1 --config samples=config.csv
```

### Key File Patterns
- **Compressed outputs**: Most analysis outputs are `.tsv.gz` files
- **BAM files**: Alignments stored in `out/aligned/{sample}.bam`
- **Variants**: `out/variants/reads/{sample}.tsv.gz` contains per-read variant calls
- **Parent assignment**: `out/parents/assigned/{sample}_assigned-parents.tsv.gz`

### Container Usage
Always specify container directive in rules:
```python
rule example:
    container: "docker://szsctt/lr_pybio:py310"
```

### Wildcard Constraints
Critical for proper rule matching with parent and sample names:
```python
wildcard_constraints:
    sample_name = '|'.join(samples.sample_name) + "|" + "|".join(samples.parent_name)
```

## Dependencies and Environment
- Primary environment defined in `mamba_env.yml`
- Core Python dependencies: biopython, pysam, numpy, scipy, pandas, tqdm
- Bioinformatics tools: minimap2, samtools (via containers)
- C3POa for consensus calling (RCA data only)

When modifying the pipeline, always consider the sequencing technology context (`seq_tech`) as it affects processing paths, especially for RCA data which requires additional consensus calling steps.