import gzip
import subprocess

import pytest
from Bio import SeqIO

from aavolve.qc_trimming_msa import (
    build_msa,
    detect_seq_format,
    parse_args,
    run_mafft_add,
    run_mafft,
    write_reference_and_reads,
)


def _write_gz_text(path: str, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_detect_seq_format_fasta_gz(tmp_path):
    path = tmp_path / "reads.fasta.gz"
    _write_gz_text(str(path), ">seq1\nACGT\n")
    assert detect_seq_format(str(path)) == "fasta"


def test_detect_seq_format_fastq_gz(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    _write_gz_text(str(path), "@seq1\nACGT\n+\n####\n")
    assert detect_seq_format(str(path)) == "fastq"


def test_detect_seq_format_empty_raises(tmp_path):
    path = tmp_path / "empty.gz"
    _write_gz_text(str(path), "\n\n")
    with pytest.raises(ValueError):
        detect_seq_format(str(path))


def test_write_reference_and_reads_limits_and_prefixes(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">refA some desc\nACGT\n>refB\nTGCA\n")

    trimmed = tmp_path / "trimmed.fastq.gz"
    _write_gz_text(
        str(trimmed),
        "@r1\nACGT\n+\n####\n@r2\nTGCA\n+\n####\n@r3\nAAAA\n+\n####\n",
    )

    combined = tmp_path / "combined.fa"
    written = write_reference_and_reads(
        reference_path=str(reference),
        reads_path=str(trimmed),
        output_fasta_path=str(combined),
        max_reads=2,
    )

    assert written == 4
    ids = [record.id for record in SeqIO.parse(str(combined), "fasta")]
    assert ids == ["ref__refA", "ref__refB", "read001__r1", "read002__r2"]


def test_run_mafft_invokes_subprocess(monkeypatch, tmp_path):
    input_fasta = tmp_path / "in.fa"
    input_fasta.write_text(">a\nACGT\n")
    output_fasta = tmp_path / "out.fa"

    called = {}

    def fake_run(cmd, check, stdout):
        called["cmd"] = cmd
        called["check"] = check
        stdout.write(">a\nACGT\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    run_mafft(input_fasta=str(input_fasta), output_fasta=str(output_fasta), threads=3)

    assert called["cmd"][:4] == ["mafft", "--auto", "--thread", "3"]
    assert called["cmd"][4] == str(input_fasta)
    assert output_fasta.read_text().startswith(">a")


def test_run_mafft_add_invokes_subprocess(monkeypatch, tmp_path):
    reads_fasta = tmp_path / "reads.fa"
    reads_fasta.write_text(">r1\nACGT\n")
    reference_alignment = tmp_path / "ref.fa"
    reference_alignment.write_text(">ref\nACGT\n")
    output_fasta = tmp_path / "out.fa"

    called = {}

    def fake_run(cmd, check, stdout):
        called["cmd"] = cmd
        called["check"] = check
        stdout.write(">ref\nACGT\n>r1\nACGT\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    run_mafft_add(
        reference_alignment_fasta=str(reference_alignment),
        reads_fasta=str(reads_fasta),
        output_fasta=str(output_fasta),
        threads=2,
    )

    assert called["cmd"][:3] == ["mafft", "--thread", "2"]
    assert "--add" in called["cmd"]
    assert output_fasta.read_text().startswith(">ref")


def test_build_msa_writes_output(monkeypatch, tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">ref\nACGT\n")

    trimmed = tmp_path / "trimmed.fasta.gz"
    _write_gz_text(str(trimmed), ">r1\nACGT\n>r2\nACGA\n")

    output = tmp_path / "out" / "msa.fa"

    called = {}

    def fake_run(cmd, check, stdout):
        called.setdefault("cmds", []).append(cmd)
        stdout.write(">ref__ref\nACGT\n>read001__r1\nACGT\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    build_msa(
        reference_path=str(reference),
        reads_path=str(trimmed),
        output_msa_path=str(output),
        max_reads=1,
        threads=2,
    )

    assert output.exists()
    assert output.read_text().startswith(">ref__ref")
    assert any("--add" in cmd for cmd in called.get("cmds", []))


def test_parse_args_accepts_reads(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">ref\nACGT\n")
    reads = tmp_path / "reads.fasta.gz"
    _write_gz_text(str(reads), ">r1\nACGT\n")
    output = tmp_path / "out.fa"

    args = parse_args(
        [
            "--reference",
            str(reference),
            "--reads",
            str(reads),
            "--output",
            str(output),
        ]
    )

    assert str(args.reference) == str(reference)
    assert str(args.reads) == str(reads)
    assert args.trimmed is None
