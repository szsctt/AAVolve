# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **BREAKING**: Upgraded from Snakemake 7 to Snakemake 8
  - Updated `snakemake-minimal` dependency from version 7 to version 8 in test environment (`tests/envs/test.yml`)
  - CI/CD workflows now use Snakemake 8.4 (`deps/snakemake/snakemake.yml`)
  - This is a major version upgrade with breaking API changes in Snakemake's internal APIs
  - Snakemake 8 introduces a new layered API architecture (`SnakemakeApi`, `DAGApi`, etc.)
  - **Users running the pipeline should install Snakemake 8.x instead of Snakemake 7.x**
  - The Snakefile itself remains compatible - no changes needed to rule syntax for end users
  - Test suite has been updated to use Snakemake 8 APIs for programmatic workflow testing
  - Removed redundant `conda:` directives from rules in favor of `container:` directives (Snakemake 8 best practice)

### Added
- Comprehensive DAG structure tests using Snakemake 8 API (`tests/snakemake/test_dag.py`)
- Anchor sequences feature for stabilizing alignments of highly similar parental sequences
- Log files for all pipeline rules

### Dependencies Updated
- `snakemake-minimal`: 7 → 8
- `pytest`: 7.4 → 9
- `quarto`: 1.4 → 1.8

### Migration Notes for Users
If you have an existing AAVolve installation with Snakemake 7:
1. Update your Snakemake installation to version 8.x:
   ```bash
   # Using conda/mamba
   mamba install -c conda-forge -c bioconda snakemake-minimal=8
   ```
2. No changes to your config files or Snakefile usage are required
3. The pipeline command-line interface remains the same

### Migration Notes for Developers
If you are developing or testing AAVolve:
1. The test environment now requires Snakemake 8.x
2. Programmatic access to Snakemake workflows now uses the new layered API:
   - `from snakemake.api import SnakemakeApi, DAGApi`
   - `from snakemake.settings.types import OutputSettings, ResourceSettings, etc.`
3. See `tests/snakemake/test_dag.py` for examples of the new API usage
4. The old Snakemake 7 API (`snakemake.snakemake()`) is deprecated in favor of the new API

### References
- [Snakemake 8 Release Notes](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html)
- [Snakemake 8 API Documentation](https://snakemake.github.io/snakemake-plugin-catalog/index.html)
