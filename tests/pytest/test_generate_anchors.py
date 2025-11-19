from Bio import SeqIO

import pytest

from aavolve.generate_anchors import parse_args, write_anchors, main
from aavolve.utils import use_open


def _read_sequences(path):
    fmt = "fasta" if path.endswith("fa") or path.endswith("fasta") else "fastq"
    with use_open(path, "rt") as handle:
        return list(SeqIO.parse(handle, fmt))


def test_write_anchors_reproducible(tmp_path):
    output = tmp_path / "anchors.fasta"
    write_anchors(5, output, seed="sample")
    records = _read_sequences(str(output))
    assert len(records) == 2
    assert records[0].id == "5prime_anchor"
    assert records[1].id == "3prime_anchor"
    assert len(records[0].seq) == 5
    assert len(records[1].seq) == 5

    second = tmp_path / "anchors_second.fasta"
    write_anchors(5, second, seed="sample")
    records_second = _read_sequences(str(second))
    assert [str(rec.seq) for rec in records] == [str(rec.seq) for rec in records_second]


def test_write_anchors_invalid_length(tmp_path):
    output = tmp_path / "anchors.fasta"
    with pytest.raises(ValueError):
        write_anchors(0, output)


def test_write_anchors_creates_parent_dirs(tmp_path):
    nested_output = tmp_path / "nested" / "anchors.fasta"
    write_anchors(3, nested_output, seed="x")
    records = _read_sequences(str(nested_output))
    assert len(records) == 2


def test_parse_args_roundtrip():
    args = parse_args(["--length", "7", "--output", "file.fa", "--seed", "seed"])
    assert args.length == 7
    assert args.output == "file.fa"
    assert args.seed == "seed"


def test_main_writes_file(tmp_path):
    out_path = tmp_path / "cli.fa"
    main(["--length", "4", "--output", str(out_path), "--seed", "seed"])
    records = _read_sequences(str(out_path))
    assert len(records) == 2
    assert all(len(rec.seq) == 4 for rec in records)
