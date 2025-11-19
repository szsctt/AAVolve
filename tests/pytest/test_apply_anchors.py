from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pytest

from aavolve.apply_anchors import (
    DEFAULT_PHRED,
    apply_anchors,
    detect_format,
    load_anchors,
    main,
    pad_quality,
    parse_args,
)
from aavolve.generate_anchors import write_anchors
from aavolve.utils import use_open


def _read_sequences(path):
    fmt = "fasta" if path.endswith("fa") or path.endswith("fasta") else "fastq"
    with use_open(path, "rt") as handle:
        return list(SeqIO.parse(handle, fmt))


def test_apply_anchors_fasta(tmp_path):
    anchor_file = tmp_path / "anchors.fasta"
    write_anchors(3, anchor_file, seed="seed")

    input_fasta = tmp_path / "input.fa"
    with use_open(input_fasta, "wt") as handle:
        handle.write(">seq1\nACGT\n")

    output_fasta = tmp_path / "anchored.fa"
    apply_anchors(str(input_fasta), str(output_fasta), str(anchor_file))

    anchors = _read_sequences(str(anchor_file))
    anchored_records = _read_sequences(str(output_fasta))
    assert len(anchored_records) == 1
    expected_prefix = str(anchors[0].seq)
    expected_suffix = str(anchors[1].seq)
    seq = str(anchored_records[0].seq)
    assert seq.startswith(expected_prefix)
    assert seq.endswith(expected_suffix)
    assert seq[len(expected_prefix):-len(expected_suffix)] == "ACGT"


def test_apply_anchors_fastq(tmp_path):
    anchor_file = tmp_path / "anchors.fasta"
    write_anchors(4, anchor_file, seed="seed")
    anchors = _read_sequences(str(anchor_file))
    five = str(anchors[0].seq)
    three = str(anchors[1].seq)

    input_fastq = tmp_path / "input.fastq"
    with use_open(input_fastq, "wt") as handle:
        handle.write("@read1\nACGT\n+\n!\"#$\n")

    output_fastq = tmp_path / "anchored.fastq"
    apply_anchors(str(input_fastq), str(output_fastq), str(anchor_file))

    anchored_records = _read_sequences(str(output_fastq))
    assert len(anchored_records) == 1
    record = anchored_records[0]
    seq = str(record.seq)
    assert seq.startswith(five)
    assert seq.endswith(three)
    assert seq[len(five):-len(three)] == "ACGT"

    qualities = record.letter_annotations["phred_quality"]
    assert len(qualities) == len(seq)
    assert qualities[:len(five)] == [DEFAULT_PHRED] * len(five)
    assert qualities[-len(three):] == [DEFAULT_PHRED] * len(three)
    assert qualities[len(five):-len(three)] == [0, 1, 2, 3]


def test_load_anchors_returns_sequences(tmp_path):
    anchor_file = tmp_path / "anchors.fa"
    with use_open(anchor_file, "wt") as handle:
        handle.write(">5prime_anchor\nacg\n>3prime_anchor\nttg\n")

    five, three = load_anchors(str(anchor_file))
    assert five == "ACG"
    assert three == "TTG"


def test_load_anchors_missing_entry(tmp_path):
    anchor_file = tmp_path / "anchors.fa"
    with use_open(anchor_file, "wt") as handle:
        handle.write(">5prime_anchor\nAAA\n")

    with pytest.raises(ValueError):
        load_anchors(str(anchor_file))


def test_detect_format(tmp_path):
    fasta = tmp_path / "input.fa"
    fastq = tmp_path / "input.fq"
    other = tmp_path / "input.txt"
    fasta.write_text(">seq\n")
    fastq.write_text("@seq\n")
    other.write_text("#")

    assert detect_format(str(fasta)) == "fasta"
    assert detect_format(str(fastq)) == "fastq"
    with pytest.raises(ValueError):
        detect_format(str(other))


def test_pad_quality_applies_defaults():
    record = SeqRecord(Seq("AGTCA"), id="r1")
    record.letter_annotations = {}
    pad_quality(record, [1, 2], 2, 1)
    assert record.letter_annotations["phred_quality"] == [DEFAULT_PHRED, DEFAULT_PHRED, 1, 2, DEFAULT_PHRED]


def test_apply_anchors_fastq_missing_phred(tmp_path, monkeypatch):
    anchor_file = tmp_path / "anchors.fasta"
    write_anchors(2, anchor_file, seed="seed")

    input_fastq = tmp_path / "input.fastq"
    input_fastq.write_text("@read1\nAC\n+\n!!\n")

    original_parse = SeqIO.parse

    def fake_parse(handle, fmt):
        if fmt == "fastq":
            record = SeqRecord(Seq("AC"), id="read1")
            record.letter_annotations = {}
            return iter([record])
        return original_parse(handle, fmt)

    monkeypatch.setattr("aavolve.apply_anchors.SeqIO.parse", fake_parse)

    output_fastq = tmp_path / "anchored.fastq"
    with pytest.raises(ValueError):
        apply_anchors(str(input_fastq), str(output_fastq), str(anchor_file))


def test_parse_args_roundtrip():
    args = parse_args(["--input", "in.fa", "--output", "out.fa", "--anchors", "anchors.fa"])
    assert args.input == "in.fa"
    assert args.output == "out.fa"
    assert args.anchors == "anchors.fa"


def test_main_creates_output(tmp_path):
    anchor_file = tmp_path / "anchors.fasta"
    write_anchors(3, anchor_file, seed="seed")
    input_fasta = tmp_path / "input.fa"
    input_fasta.write_text(">seq\nACGT\n")
    output_fasta = tmp_path / "out.fa"

    main([
        "--input",
        str(input_fasta),
        "--output",
        str(output_fasta),
        "--anchors",
        str(anchor_file),
    ])

    records = _read_sequences(str(output_fasta))
    assert len(records) == 1
