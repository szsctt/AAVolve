import pytest


pysam = pytest.importorskip("pysam")


def _write_fasta(path, records):
    with open(path, "wt") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def _write_parent_bam(path, *, ref_name="REF", ref_len=10, alignments):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": ref_name, "LN": ref_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for qname, start, cigar in alignments:
            a = pysam.AlignedSegment()
            a.query_name = qname
            qlen = int(cigar[:-1]) if cigar.endswith("M") else ref_len
            a.query_sequence = "A" * qlen
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
            bam.write(a)
    pysam.index(str(path))


def test_warn_parents_dropped_empty_when_all_ok(tmp_path):
    from aavolve.warn_parents_dropped import main

    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"
    bam = tmp_path / "parents.bam"
    out = tmp_path / "warn.txt"

    _write_fasta(parents, [("P1", "A" * 10), ("P2", "A" * 10)])
    _write_fasta(reference, [("REF", "A" * 10)])
    _write_parent_bam(bam, alignments=[("P1", 0, "10M"), ("P2", 0, "10M")])

    rc = main(
        [
            "--parents-fasta",
            str(parents),
            "--bam",
            str(bam),
            "--reference-fasta",
            str(reference),
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    assert out.read_text() == ""


def test_warn_parents_dropped_reports_partial_alignment(tmp_path):
    from aavolve.warn_parents_dropped import main

    parents = tmp_path / "parents.fa"
    reference = tmp_path / "ref.fa"
    bam = tmp_path / "parents.bam"
    out = tmp_path / "warn.txt"

    _write_fasta(parents, [("P1", "A" * 10), ("P2", "A" * 10)])
    _write_fasta(reference, [("REF", "A" * 10)])
    _write_parent_bam(bam, alignments=[("P1", 0, "10M"), ("P2", 0, "9M")])  # P2 not end-to-end

    rc = main(
        [
            "--parents-fasta",
            str(parents),
            "--bam",
            str(bam),
            "--reference-fasta",
            str(reference),
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    text = out.read_text()
    assert "dropped parents" in text.lower()
    assert "P2" in text

