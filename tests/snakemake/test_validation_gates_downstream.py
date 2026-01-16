from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest


snakemake = pytest.importorskip("snakemake")


@pytest.fixture(scope="session")
def snakefile() -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    p = repo_root / "Snakefile"
    if not p.exists():
        raise FileNotFoundError(f"Could not find Snakefile at expected location: {p}")
    return p


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    lines = []
    for name, seq in records:
        lines.append(f">{name}\n{seq}\n")
    path.write_text("".join(lines))


def test_failed_validate_inputs_blocks_downstream_jobs(tmp_path: Path, snakefile: Path) -> None:
    """
    If validate_inputs fails, the workflow should not start any downstream work.

    This test runs a small workflow in a temp workdir with an intentionally-invalid
    parents FASTA and asserts that a representative downstream output is not created.
    """
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fastq"
    samples_csv = tmp_path / "samples.csv"

    # Trigger a validate_inputs failure: duplicate parent IDs.
    _write_fasta(parents, [("P1", "ACGT"), ("P1", "ACGA")])
    _write_fasta(reference, [("REF", "ACGT")])
    reads.write_text("@r1\nACGT\n+\n####\n")

    samples_csv.write_text(
        "sample_name,parent_name,reference_name,seq_tech,read_file,parent_file,reference_file\n"
        f"s1,p1,ref,np,{reads},{parents},{reference}\n"
    )

    env = os.environ.copy()
    # Ensure aavolve is importable for the Python invocations inside rules.
    env["PYTHONPATH"] = os.pathsep.join([str(snakefile.parent), env.get("PYTHONPATH", "")]).strip(
        os.pathsep
    )
    # Avoid writing to a non-writable $HOME cache in sandboxed/CI environments.
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(tmp_path),
            "--cores",
            "1",
            "--config",
            f"samples={samples_csv}",
        ],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode != 0
    # The validation error should be visible.
    assert "input validation failed" in (result.stderr + result.stdout).lower()

    # Downstream work should not start, so representative outputs should not exist.
    assert not (tmp_path / "out" / "aligned" / "s1.bam").exists()
    assert not (tmp_path / "out" / "qc" / "s1_read-counts.tsv").exists()
    assert not (tmp_path / "out" / "qc" / "input-checks" / "_all.ok").exists()
