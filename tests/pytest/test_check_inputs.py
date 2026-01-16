import gzip
import subprocess
import sys
from pathlib import Path
import os

import pytest

from aavolve.check_inputs import validate_inputs


def _write_fasta(path, records):
    lines = []
    for name, seq in records:
        lines.append(f">{name}\n{seq}\n")
    path.write_text("".join(lines))


def _write_fasta_gz(path, records):
    lines = []
    for name, seq in records:
        lines.append(f">{name}\n{seq}\n")
    with gzip.open(path, "wt") as handle:
        handle.write("".join(lines))


def test_validate_inputs_happy_path(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT"), ("P2", "ACGA")])
    _write_fasta(reference, [("REF", "ACGT")])

    validate_inputs(str(parents), str(reference))


def test_validate_inputs_reference_must_be_single_record(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT")])
    _write_fasta(reference, [("R1", "ACGT"), ("R2", "ACGT")])

    with pytest.raises(ValueError, match="exactly one sequence"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_parent_ids_unique(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT"), ("P1", "ACGT")])
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match="Duplicate sequence name"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_no_duplicate_parent_sequences(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT"), ("P2", "ACGT")])
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match="Duplicate parent sequences"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_reference_parent_same_name_same_sequence(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("REF", "ACGT"), ("P2", "ACGA")])
    _write_fasta(reference, [("REF", "ACGA")])

    with pytest.raises(ValueError, match="sequences differ"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_accepts_gz(tmp_path):
    parents = tmp_path / "parents.fa.gz"
    reference = tmp_path / "ref.fa.gz"

    _write_fasta_gz(parents, [("P1", "ACGT"), ("P2", "ACGA")])
    _write_fasta_gz(reference, [("REF", "ACGT")])

    validate_inputs(str(parents), str(reference))


def test_validate_inputs_missing_parents_file(tmp_path):
    parents = tmp_path / "nope.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match="Parents FASTA not found"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_missing_reference_file(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "nope.fa"

    _write_fasta(parents, [("P1", "ACGT")])

    with pytest.raises(ValueError, match="Reference FASTA not found"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_parents_empty_file(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    parents.write_text("")
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match="Parents FASTA is empty"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_reference_empty_file(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT")])
    reference.write_text("")

    with pytest.raises(ValueError, match="Reference FASTA is empty"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_rejects_bad_parent_id_whitespace(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P 1", "ACGT")])
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match="whitespace or commas"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_rejects_bad_parent_id_comma(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P,1", "ACGT")])
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match="whitespace or commas"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_rejects_empty_parent_sequence(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "")])
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match=r"Parent sequence P1 is empty"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_rejects_invalid_parent_bases(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGTN"), ("P2", "ACGTX")])
    _write_fasta(reference, [("REF", "ACGT")])

    with pytest.raises(ValueError, match=r"Parent sequence P2 has invalid base"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_rejects_empty_reference_sequence(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT")])
    _write_fasta(reference, [("REF", "")])

    with pytest.raises(ValueError, match="Reference sequence is empty"):
        validate_inputs(str(parents), str(reference))


def test_validate_inputs_rejects_invalid_reference_bases(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"

    _write_fasta(parents, [("P1", "ACGT")])
    _write_fasta(reference, [("REF", "ACGTX")])

    with pytest.raises(ValueError, match=r"Reference sequence has invalid base"):
        validate_inputs(str(parents), str(reference))


def test_cli_writes_ok_file_on_success(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"
    output = tmp_path / "out" / "ok.txt"

    _write_fasta(parents, [("P1", "ACGT"), ("P2", "ACGA")])
    _write_fasta(reference, [("REF", "ACGT")])

    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(Path(__file__).resolve().parents[2]), env.get("PYTHONPATH", "")]
    ).strip(os.pathsep)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "aavolve.check_inputs",
            "--parents",
            str(parents),
            "--reference",
            str(reference),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stderr
    assert output.read_text() == "ok\n"


def test_cli_exits_2_and_prints_error_on_failure(tmp_path):
    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"
    output = tmp_path / "out" / "ok.txt"

    _write_fasta(parents, [("P1", "ACGT"), ("P1", "ACGT")])
    _write_fasta(reference, [("REF", "ACGT")])

    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(Path(__file__).resolve().parents[2]), env.get("PYTHONPATH", "")]
    ).strip(os.pathsep)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "aavolve.check_inputs",
            "--parents",
            str(parents),
            "--reference",
            str(reference),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 2
    assert "ERROR:" in result.stderr
    assert not output.exists()
