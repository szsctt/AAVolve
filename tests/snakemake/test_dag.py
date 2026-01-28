from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from snakemake.api import SnakemakeApi, DAGApi
from snakemake.settings.types import (
    OutputSettings,
    ResourceSettings,
    StorageSettings,
    ConfigSettings,
    DAGSettings,
)


from aavolve.get_samples import get_samples

# Skip tests if snakemake isn't available
snakemake = pytest.importorskip("snakemake")

@pytest.fixture(scope="session")
def snakefile() -> Path:
    """
    Resolve the Snakefile path: <repo root>/Snakefile.
    """
    # Assume repository root has Snakefile
    repo_root = Path(__file__).resolve().parents[2]
    p = repo_root / "Snakefile"
    if not p.exists():
        raise FileNotFoundError(f"Could not find Snakefile at expected location: {p}")
    return p

@pytest.fixture
def anchors_cc_config() -> Path:
    """
    Provide path to a minimal config CSV for testing anchor-related rules.
    """
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "tests" / "data" / "config" / "test_workflow_cc.csv"
    if not config_path.exists():
        raise FileNotFoundError(f"Could not find test config at expected location: {config_path}")
    return config_path

@pytest.fixture
def samples_config() -> Path:
    """
    Provide path to test samples config with np-cc samples.
    """
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "tests" / "data" / "config" / "test_samples.csv"
    if not config_path.exists():
        raise FileNotFoundError(f"Could not find test config at expected location: {config_path}")
    return config_path

@pytest.fixture
def anchors_config() -> Path:
    """
    Provide path to test config with anchors enabled.
    """
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "tests" / "data" / "config" / "test_anchors.csv"
    if not config_path.exists():
        raise FileNotFoundError(f"Could not find test config at expected location: {config_path}")
    return config_path

@pytest.fixture
def np_only_config() -> Path:
    """
    Provide path to test config with only non-np-cc (np) samples.
    """
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "tests" / "data" / "config" / "test_np_only.csv"
    if not config_path.exists():
        raise FileNotFoundError(f"Could not find test config at expected location: {config_path}")
    return config_path


@pytest.fixture
def npcc_dedup_minreps_config() -> Path:
    """Test config where np-cc samples share reads/splint but differ in min_reps."""
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "tests" / "data" / "config" / "test_npcc_dedup_minreps.csv"
    if not config_path.exists():
        raise FileNotFoundError(f"Could not find test config at expected location: {config_path}")
    return config_path

def build_dag(
    snakefile: Path,
    targets: list[str],
    config: dict,
    cores: int = 1,
) -> DAGApi:
    """
    Construct a DAGApi instance using the Snakemake 8+ layered API.
    Returns the DAGApi object for inspection (no execution).
    
    Note: workdir is set to the repo root so that relative paths in config resolve correctly.
    """
    # Use repo root as workdir so relative paths in test configs work
    repo_root = snakefile.parent
    
    with SnakemakeApi(OutputSettings(verbose=False)) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            resource_settings=ResourceSettings(cores=cores),
            config_settings=ConfigSettings(config=config),
            storage_settings=StorageSettings(),
            snakefile=snakefile,
            workdir=repo_root,
        )
        
        dag_api = workflow_api.dag(
            dag_settings=DAGSettings(targets=targets)
        )
        
        return dag_api

class DAGContext:
    """
    Context manager for building and inspecting a Snakemake DAG with access to expanded jobs.
    
    Usage:
        with build_dag_with_jobs(snakefile, config) as dag_context:
            workflow = dag_context.workflow
            jobs = dag_context.jobs
            # inspect jobs with expanded wildcards
    """
    def __init__(self, snakefile: Path, config: dict, targets: list[str] = None, cores: int = 1):
        self.snakefile = snakefile
        self.config = config
        self.targets = targets or ["all"]
        self.cores = cores
        self.repo_root = snakefile.parent
        self._snakemake_api = None
        self._workflow_api = None
        self._dag_api = None
        self.workflow = None
        self.jobs = None
    
    def __enter__(self):
        self._snakemake_api = SnakemakeApi(OutputSettings(verbose=False))
        self._snakemake_api.__enter__()
        
        self._workflow_api = self._snakemake_api.workflow(
            snakefile=self.snakefile,
            workdir=self.repo_root,
            resource_settings=ResourceSettings(cores=self.cores),
            config_settings=ConfigSettings(config=self.config),
            storage_settings=StorageSettings(),
        )
        
        self._dag_api = self._workflow_api.dag(
            dag_settings=DAGSettings(targets=self.targets)
        )
        
        # Trigger DAG build by calling summary() (which populates workflow.dag)
        self._dag_api.summary()
        
        # Access the workflow and its jobs
        self.workflow = self._dag_api.workflow_api._workflow
        if self.workflow.dag is not None:
            self.jobs = list(self.workflow.dag.jobs)
        else:
            self.jobs = []
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._snakemake_api:
            self._snakemake_api.__exit__(exc_type, exc_val, exc_tb)
        return False

def test_rulegraph_contains_expected_rules(tmp_path, snakefile, anchors_cc_config):
    """Test that the DAG contains expected rules when built with the test config."""
    dag_api = build_dag(
        snakefile=snakefile,
        targets=["all"],
        config={"samples": str(anchors_cc_config)},
    )

    # Access the workflow to check rules
    workflow = dag_api.workflow_api._workflow
    rule_names = {r.name for r in workflow.rules}
    
    # Assert presence of key rules for the pipeline
    assert "all" in rule_names, "Missing 'all' rule"
    assert "consensus" in rule_names, "Missing 'consensus' rule"
    assert "filter_consensus" in rule_names, "Missing 'filter_consensus' rule"
    assert "align" in rule_names, "Missing 'align' rule"
    assert "extract_variants_reads" in rule_names, "Missing 'extract_variants_reads' rule"
    assert "assign_parents" in rule_names, "Missing 'assign_parents' rule"
    
    # Check that consensus rule has expected output pattern
    consensus_rule = workflow.get_rule("consensus")
    assert "out/c3poa/" in str(consensus_rule.output), "Consensus rule output doesn't match expected pattern"

def test_job_dependency_expansion(tmp_path, snakefile, np_only_config):
    """Test that jobs are expanded correctly with actual file paths for each sample."""
    
    # Load samples to verify expansion
    samples_df = get_samples({"samples": str(np_only_config)})
    assert len(samples_df) > 0, "Test config should have at least one sample"
    sample_name = samples_df.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(np_only_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        assert len(jobs_list) > 0, "DAG should contain jobs"
        
        # Find align job for our sample
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{sample_name}'"
        align_job = align_jobs[0]
        
        # Check actual expanded output paths
        align_outputs = list(align_job.output)
        assert f"out/aligned/{sample_name}.bam" in align_outputs, \
            f"Align job should output 'out/aligned/{sample_name}.bam', got: {align_outputs}"
        assert f"out/aligned/{sample_name}.bam.bai" in align_outputs, \
            f"Align job should output 'out/aligned/{sample_name}.bam.bai', got: {align_outputs}"
        
        # Find extract_variants_reads job
        extract_jobs = [j for j in jobs_list if j.rule.name == 'extract_variants_reads' and dict(j.wildcards).get('sample') == sample_name]
        assert len(extract_jobs) > 0, f"Should have extract_variants_reads job for sample '{sample_name}'"
        extract_outputs = list(extract_jobs[0].output)
        assert f"out/variants/reads/{sample_name}.tsv.gz" in extract_outputs, \
            f"Extract variants job should output 'out/variants/reads/{sample_name}.tsv.gz', got: {extract_outputs}"
        
        # Find assign_parents job
        assign_jobs = [j for j in jobs_list if j.rule.name == 'assign_parents' and dict(j.wildcards).get('sample') == sample_name]
        assert len(assign_jobs) > 0, f"Should have assign_parents job for sample '{sample_name}'"
        assign_outputs = list(assign_jobs[0].output)
        assert f"out/parents/assigned/{sample_name}_assigned-parents.tsv.gz" in assign_outputs, \
            f"Assign parents job should output 'out/parents/assigned/{sample_name}_assigned-parents.tsv.gz', got: {assign_outputs}"


def test_group_manifest_depends_on_samples_csv(tmp_path, snakefile, np_only_config):
    with DAGContext(snakefile, {"samples": str(np_only_config)}) as dag_ctx:
        manifest_jobs = [j for j in dag_ctx.jobs if j.rule.name == "group_manifest"]
        assert manifest_jobs, "Expected at least one group_manifest job in the DAG"
        for job in manifest_jobs:
            assert str(np_only_config) in {str(p) for p in job.input}, (
                "group_manifest should depend on the samples CSV so it is regenerated when samples change"
            )



def test_consensus_rules_present(tmp_path, snakefile, samples_config):
    """Test that consensus-related rules are present for np-cc samples."""
    # Verify that config has np-cc samples
    samples_df = get_samples({"samples": str(samples_config)})
    npcc_samples = samples_df[samples_df.seq_tech == 'np-cc']
    assert len(npcc_samples) > 0, "samples_config should have np-cc samples for this test"
    npcc_sample_name = npcc_samples.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(samples_config)}) as dag_ctx:
        workflow = dag_ctx.workflow
        rule_names = {r.name for r in workflow.rules}
        
        # Verify consensus-specific rules are present
        assert "consensus" in rule_names, "Missing 'consensus' rule"
        assert "filter_consensus" in rule_names, "Missing 'filter_consensus' rule"
        
        jobs_list = dag_ctx.jobs
        
        # Find consensus job for np-cc sample and check actual expanded paths
        consensus_jobs = [j for j in jobs_list if j.rule.name == 'consensus' and dict(j.wildcards).get('sample') == npcc_sample_name]
        assert len(consensus_jobs) > 0, f"Should have consensus job for np-cc sample '{npcc_sample_name}'"
        consensus_outputs = list(consensus_jobs[0].output)
        assert f"out/c3poa/{npcc_sample_name}/splint/R2C2_Consensus.fasta.gz" in consensus_outputs, \
            f"Consensus job should output 'out/c3poa/{npcc_sample_name}/splint/R2C2_Consensus.fasta.gz', got: {consensus_outputs}"
        
        # Find filter_consensus job and check actual expanded paths
        filter_jobs = [j for j in jobs_list if j.rule.name == 'filter_consensus' and dict(j.wildcards).get('sample') == npcc_sample_name]
        assert len(filter_jobs) > 0, f"Should have filter_consensus job for np-cc sample '{npcc_sample_name}'"
        filter_outputs = list(filter_jobs[0].output)
        assert f"out/c3poa_filt/{npcc_sample_name}.fasta.gz" in filter_outputs, \
            f"Filter consensus job should output 'out/c3poa_filt/{npcc_sample_name}.fasta.gz', got: {filter_outputs}"


def test_npcc_consensus_deduplicated_by_inputs(tmp_path, snakefile, samples_config):
    """Consensus should run once per unique (read_file,splint_file) pair, even if multiple samples share inputs."""
    samples_df = get_samples({"samples": str(samples_config)})
    npcc_samples = samples_df[samples_df.seq_tech == "np-cc"]
    assert len(npcc_samples) > 0, "samples_config should have np-cc samples for this test"
    unique_ids = sorted(set(npcc_samples["npcc_consensus_id"].astype(str)))

    with DAGContext(snakefile, {"samples": str(samples_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        dedup_jobs = [j for j in jobs_list if j.rule.name == "consensus_by_input"]
        assert len(dedup_jobs) == len(unique_ids), (
            f"Expected {len(unique_ids)} consensus_by_input jobs (one per unique input), got {len(dedup_jobs)}"
        )
        job_ids = sorted({dict(j.wildcards).get("cc_id") for j in dedup_jobs})
        assert job_ids == unique_ids


def test_npcc_consensus_dedup_when_minreps_differs(tmp_path, snakefile, npcc_dedup_minreps_config):
    """Different min_reps should NOT cause consensus_by_input to rerun; only filter_consensus differs per sample."""
    samples_df = get_samples({"samples": str(npcc_dedup_minreps_config)})
    npcc_samples = samples_df[samples_df.seq_tech == "np-cc"]
    assert len(npcc_samples) == 2
    assert len(set(npcc_samples["npcc_consensus_id"].astype(str))) == 1, "Inputs are shared; expected 1 consensus id"

    with DAGContext(snakefile, {"samples": str(npcc_dedup_minreps_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        consensus_jobs = [j for j in jobs_list if j.rule.name == "consensus_by_input"]
        assert len(consensus_jobs) == 1, "Consensus should run once per shared (read_file,splint_file) pair"

        filter_jobs = [j for j in jobs_list if j.rule.name == "filter_consensus"]
        assert len(filter_jobs) == 2, "Filtering should run per sample because min_reps differs"
        filter_outputs = sorted(str(p) for j in filter_jobs for p in j.output)
        assert filter_outputs == sorted(
            [f"out/c3poa_filt/{name}.fasta.gz" for name in npcc_samples["sample_name"].tolist()]
        )
def test_consensus_not_performed_for_non_npcc(tmp_path, snakefile, np_only_config):
    """Test that consensus jobs are not created for non-np-cc samples."""
    dag_api = build_dag(
        snakefile=snakefile,
        targets=["all"],
        config={"samples": str(np_only_config)},
    )
    
    # Access the workflow
    workflow = dag_api.workflow_api._workflow
    
    # Verify that config has only non-np-cc samples
    samples_df = get_samples({"samples": str(np_only_config)})
    assert len(samples_df) > 0, "Test config should have at least one sample"
    npcc_samples = samples_df[samples_df.seq_tech == 'np-cc']
    assert len(npcc_samples) == 0, "Test config should have NO np-cc samples"
    
    # Check that no consensus or filter_consensus jobs are in the DAG
    if workflow.dag:
        job_rules = {job.rule.name for job in workflow.dag.jobs}
        assert "consensus" not in job_rules, \
            "consensus job should NOT be in DAG for non-np-cc samples"
        assert "filter_consensus" not in job_rules, \
            "filter_consensus job should NOT be in DAG for non-np-cc samples"
        
        # Verify that align jobs exist (pipeline should still work)
        assert "align" in job_rules, \
            "align job should be in DAG for non-np-cc samples"


def test_trimming_when_enabled(tmp_path, snakefile, anchors_cc_config):
    """Test that trim_reads rule is present when trimming is enabled."""
    # Verify that config has trimming enabled
    samples_df = get_samples({"samples": str(anchors_cc_config)})
    trim_samples = samples_df[samples_df.trim == True]
    assert len(trim_samples) > 0, "Test config should have samples with trimming enabled"
    trim_sample_name = trim_samples.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(anchors_cc_config)}) as dag_ctx:
        workflow = dag_ctx.workflow
        rule_names = {r.name for r in workflow.rules}
        
        # Verify trim_reads rule is present
        assert "trim_reads" in rule_names, "Missing 'trim_reads' rule when trimming is enabled"
        
        jobs_list = dag_ctx.jobs
        
        # Find trim_reads job for sample with trimming enabled and check actual expanded paths
        trim_jobs = [j for j in jobs_list if j.rule.name == 'trim_reads' and dict(j.wildcards).get('sample') == trim_sample_name]
        assert len(trim_jobs) > 0, f"Should have trim_reads job for sample '{trim_sample_name}' with trimming enabled"
        trim_outputs = list(trim_jobs[0].output)
        assert f"out/trimmed/{trim_sample_name}.trimmed.gz" in trim_outputs, \
            f"Trim reads job should output 'out/trimmed/{trim_sample_name}.trimmed.gz', got: {trim_outputs}"


def test_trimming_not_required_when_disabled(tmp_path, snakefile, samples_config):
    """Test that trim_reads jobs are not created when trimming is disabled."""
    # Verify that config has samples without trimming
    samples_df = get_samples({"samples": str(samples_config)})
    no_trim_samples = samples_df[samples_df.trim == False]
    assert len(no_trim_samples) > 0, "Test config should have samples with trimming disabled"
    no_trim_sample_name = no_trim_samples.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(samples_config)}) as dag_ctx:
        workflow = dag_ctx.workflow
        rule_names = {r.name for r in workflow.rules}
        
        # The trim_reads rule is always defined
        assert "trim_reads" in rule_names, "trim_reads rule should be defined"
        
        jobs_list = dag_ctx.jobs
        
        # Verify that NO trim_reads job exists for samples with trimming disabled
        trim_jobs = [j for j in jobs_list if j.rule.name == 'trim_reads' and dict(j.wildcards).get('sample') == no_trim_sample_name]
        assert len(trim_jobs) == 0, \
            f"Should NOT have trim_reads job for sample '{no_trim_sample_name}' with trimming disabled, but found: {[list(j.output) for j in trim_jobs]}"
        
        # Verify that align jobs still exist (pipeline should still work without trimming)
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == no_trim_sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{no_trim_sample_name}' even without trimming"

def test_anchor_rules_when_enabled(tmp_path, snakefile, anchors_config):
    """Test that anchor-related rules are present when anchors are enabled."""
    # Verify that config has anchors enabled
    samples_df = get_samples({"samples": str(anchors_config)})
    # Check that anchors column has non-empty, non-zero values
    anchor_samples = samples_df[samples_df.anchors.notna() & (samples_df.anchors != '') & (samples_df.anchors != 0)]
    assert len(anchor_samples) > 0, "Test config should have samples with anchors enabled"
    anchor_sample_name = anchor_samples.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(anchors_config)}) as dag_ctx:
        workflow = dag_ctx.workflow
        rule_names = {r.name for r in workflow.rules}
        
        # Verify anchor-related rules are present
        assert "generate_anchor_sequences" in rule_names, "Missing 'generate_anchor_sequences' rule"
        assert "anchor_reads" in rule_names, "Missing 'anchor_reads' rule"
        assert "anchor_reference" in rule_names, "Missing 'anchor_reference' rule"
        
        jobs_list = dag_ctx.jobs
        
        # Find generate_anchor_sequences job and check actual expanded paths
        gen_anchor_jobs = [j for j in jobs_list if j.rule.name == 'generate_anchor_sequences' and dict(j.wildcards).get('sample') == anchor_sample_name]
        assert len(gen_anchor_jobs) > 0, f"Should have generate_anchor_sequences job for sample '{anchor_sample_name}' with anchors enabled"
        gen_anchor_outputs = list(gen_anchor_jobs[0].output)
        assert f"out/anchors/{anchor_sample_name}.fasta" in gen_anchor_outputs, \
            f"Generate anchor sequences job should output 'out/anchors/{anchor_sample_name}.fasta', got: {gen_anchor_outputs}"
        
        # Find anchor_reads job and check actual expanded paths
        anchor_reads_jobs = [j for j in jobs_list if j.rule.name == 'anchor_reads' and dict(j.wildcards).get('sample') == anchor_sample_name]
        assert len(anchor_reads_jobs) > 0, f"Should have anchor_reads job for sample '{anchor_sample_name}' with anchors enabled"
        anchor_reads_outputs = list(anchor_reads_jobs[0].output)
        # Check that output path starts with out/anchors/reads/{sample}
        assert any(out.startswith(f"out/anchors/reads/{anchor_sample_name}") for out in anchor_reads_outputs), \
            f"Anchor reads job should output files under 'out/anchors/reads/{anchor_sample_name}', got: {anchor_reads_outputs}"
        
        # Find anchor_reference job and check actual expanded paths
        anchor_ref_jobs = [j for j in jobs_list if j.rule.name == 'anchor_reference' and dict(j.wildcards).get('sample') == anchor_sample_name]
        assert len(anchor_ref_jobs) > 0, f"Should have anchor_reference job for sample '{anchor_sample_name}' with anchors enabled"
        anchor_ref_outputs = list(anchor_ref_jobs[0].output)
        assert f"out/anchors/references/{anchor_sample_name}.fasta" in anchor_ref_outputs, \
            f"Anchor reference job should output 'out/anchors/references/{anchor_sample_name}.fasta', got: {anchor_ref_outputs}"

def test_anchor_rules_not_required_when_disabled(tmp_path, snakefile, samples_config):
    """Test that anchor jobs are not created when anchors are disabled."""
    # Verify that config has samples without anchors (or with anchors=0 or empty)
    samples_df = get_samples({"samples": str(samples_config)})
    # Most samples in test_samples.csv don't have anchors or have empty anchors column
    no_anchor_samples = samples_df[samples_df.anchors.isna() | (samples_df.anchors == '') | (samples_df.anchors == 0)]
    assert len(no_anchor_samples) > 0, "Test config should have samples without anchors"
    no_anchor_sample_name = no_anchor_samples.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(samples_config)}) as dag_ctx:
        workflow = dag_ctx.workflow
        rule_names = {r.name for r in workflow.rules}
        
        # The anchor rules are always defined
        assert "generate_anchor_sequences" in rule_names, "generate_anchor_sequences rule should be defined"
        assert "anchor_reads" in rule_names, "anchor_reads rule should be defined"
        assert "anchor_reference" in rule_names, "anchor_reference rule should be defined"
        
        jobs_list = dag_ctx.jobs
        
        # Verify that NO anchor-related jobs exist for samples with anchors disabled
        gen_anchor_jobs = [j for j in jobs_list if j.rule.name == 'generate_anchor_sequences' and dict(j.wildcards).get('sample') == no_anchor_sample_name]
        assert len(gen_anchor_jobs) == 0, \
            f"Should NOT have generate_anchor_sequences job for sample '{no_anchor_sample_name}' without anchors, but found: {[list(j.output) for j in gen_anchor_jobs]}"
        
        anchor_reads_jobs = [j for j in jobs_list if j.rule.name == 'anchor_reads' and dict(j.wildcards).get('sample') == no_anchor_sample_name]
        assert len(anchor_reads_jobs) == 0, \
            f"Should NOT have anchor_reads job for sample '{no_anchor_sample_name}' without anchors, but found: {[list(j.output) for j in anchor_reads_jobs]}"
        
        anchor_ref_jobs = [j for j in jobs_list if j.rule.name == 'anchor_reference' and dict(j.wildcards).get('sample') == no_anchor_sample_name]
        assert len(anchor_ref_jobs) == 0, \
            f"Should NOT have anchor_reference job for sample '{no_anchor_sample_name}' without anchors, but found: {[list(j.output) for j in anchor_ref_jobs]}"
        
        # Verify that align jobs still exist (pipeline should still work without anchors)
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == no_anchor_sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{no_anchor_sample_name}' even without anchors"

def test_align_inputs_with_trimming_and_anchors(tmp_path, snakefile, anchors_cc_config):
    """Test that align job uses trimmed and anchored reads when both are enabled."""
    # Get a sample with both trimming and anchors enabled
    samples_df = get_samples({"samples": str(anchors_cc_config)})
    samples_with_both = samples_df[
        (samples_df.trim == True) & 
        (samples_df.anchors.notna()) & 
        (samples_df.anchors != '') & 
        (samples_df.anchors != 0)
    ]
    assert len(samples_with_both) > 0, "Test config should have samples with both trimming and anchors enabled"
    sample_name = samples_with_both.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(anchors_cc_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        
        # Find align job
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{sample_name}'"
        align_inputs = list(align_jobs[0].input)
        
        # Align should use anchored reads (which come after trimming)
        # Check for anchored reads path (extension may vary)
        assert any(f"out/anchors/reads/{sample_name}" in str(inp) for inp in align_inputs), \
            f"Align job should use anchored reads under 'out/anchors/reads/{sample_name}' when both trimming and anchors are enabled, got: {align_inputs}"
        
        # Align should use anchored reference
        expected_ref_input = f"out/anchors/references/{sample_name}.fasta"
        assert expected_ref_input in align_inputs, \
            f"Align job should use anchored reference '{expected_ref_input}' when anchors are enabled, got: {align_inputs}"

def test_align_inputs_with_trimming_no_anchors(tmp_path, snakefile, samples_config):
    """Test that align job uses trimmed reads (not anchored) when only trimming is enabled."""
    # Get a sample with trimming enabled but no anchors
    samples_df = get_samples({"samples": str(samples_config)})
    samples_trim_no_anchor = samples_df[
        (samples_df.trim == True) & 
        (samples_df.anchors.isna() | (samples_df.anchors == '') | (samples_df.anchors == 0))
    ]
    assert len(samples_trim_no_anchor) > 0, "Test config should have samples with trimming but no anchors"
    sample_name = samples_trim_no_anchor.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(samples_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        
        # Find align job
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{sample_name}'"
        align_inputs = list(align_jobs[0].input)
        
        # Align should use trimmed reads
        expected_reads_input = f"out/trimmed/{sample_name}.trimmed.gz"
        assert expected_reads_input in align_inputs, \
            f"Align job should use trimmed reads '{expected_reads_input}' when trimming is enabled but no anchors, got: {align_inputs}"
        
        # Should NOT use anchored reads or reference
        assert not any("anchors/reads" in str(inp) for inp in align_inputs), \
            f"Align job should NOT use anchored reads when anchors disabled, got: {align_inputs}"
        assert not any("anchors/references" in str(inp) for inp in align_inputs), \
            f"Align job should NOT use anchored reference when anchors disabled, got: {align_inputs}"

def test_align_inputs_with_anchors_no_trimming(tmp_path, snakefile, anchors_config):
    """Test that align job uses anchored reads (from raw reads) when only anchors are enabled."""
    # Get a sample with anchors enabled but no trimming
    samples_df = get_samples({"samples": str(anchors_config)})
    samples_anchor_no_trim = samples_df[
        (samples_df.trim == False) & 
        (samples_df.anchors.notna()) & 
        (samples_df.anchors != '') & 
        (samples_df.anchors != 0)
    ]
    assert len(samples_anchor_no_trim) > 0, "Test config should have samples with anchors but no trimming"
    sample_name = samples_anchor_no_trim.iloc[0]['sample_name']
    
    with DAGContext(snakefile, {"samples": str(anchors_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        
        # Find align job
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{sample_name}'"
        align_inputs = list(align_jobs[0].input)
        
        # Align should use anchored reads
        assert any(f"out/anchors/reads/{sample_name}" in str(inp) for inp in align_inputs), \
            f"Align job should use anchored reads under 'out/anchors/reads/{sample_name}' when anchors are enabled, got: {align_inputs}"
        
        # Align should use anchored reference
        expected_ref_input = f"out/anchors/references/{sample_name}.fasta"
        assert expected_ref_input in align_inputs, \
            f"Align job should use anchored reference '{expected_ref_input}' when anchors are enabled, got: {align_inputs}"
        
        # Should NOT use trimmed reads
        assert not any("trimmed" in str(inp) for inp in align_inputs), \
            f"Align job should NOT use trimmed reads when trimming disabled, got: {align_inputs}"

def test_align_inputs_no_trimming_no_anchors(tmp_path, snakefile, np_only_config):
    """Test that align job uses raw reads when neither trimming nor anchors are enabled."""
    # Get a sample with no trimming and no anchors
    samples_df = get_samples({"samples": str(np_only_config)})
    sample_name = samples_df.iloc[0]['sample_name']
    
    # Verify no trimming and no anchors
    assert not samples_df.iloc[0]['trim'], "Sample should have trimming disabled"
    assert not samples_df.iloc[0]['anchors'] or pd.isna(samples_df.iloc[0]['anchors']), \
        "Sample should have anchors disabled"
    
    with DAGContext(snakefile, {"samples": str(np_only_config)}) as dag_ctx:
        jobs_list = dag_ctx.jobs
        
        # Find align job
        align_jobs = [j for j in jobs_list if j.rule.name == 'align' and dict(j.wildcards).get('sample') == sample_name]
        assert len(align_jobs) > 0, f"Should have align job for sample '{sample_name}'"
        align_inputs = list(align_jobs[0].input)
        
        # Should NOT use trimmed or anchored reads
        assert not any("trimmed" in str(inp) for inp in align_inputs), \
            f"Align job should NOT use trimmed reads when trimming disabled, got: {align_inputs}"
        assert not any("anchors/reads" in str(inp) for inp in align_inputs), \
            f"Align job should NOT use anchored reads when anchors disabled, got: {align_inputs}"
        assert not any("anchors/references" in str(inp) for inp in align_inputs), \
            f"Align job should NOT use anchored reference when anchors disabled, got: {align_inputs}"
        
        # Should use raw read file (either original or filtered consensus depending on seq_tech)
        # For np samples, should use the read_file directly or filtered consensus for np-cc
        seq_tech = samples_df.iloc[0]['seq_tech']
        if seq_tech == 'np-cc':
            expected_input_pattern = f"out/c3poa_filt/{sample_name}.fasta.gz"
            assert expected_input_pattern in align_inputs, \
                f"Align job for np-cc should use filtered consensus '{expected_input_pattern}', got: {align_inputs}"
        else:
            # For non-np-cc, should use original read file (from data/ directory typically)
            assert any("data/" in str(inp) or inp.endswith(".fastq.gz") or inp.endswith(".fasta.gz") 
                      for inp in align_inputs), \
                f"Align job should use raw read file from data directory, got: {align_inputs}"
