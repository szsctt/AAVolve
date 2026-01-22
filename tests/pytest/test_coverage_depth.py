import gzip

import pytest


pysam = pytest.importorskip("pysam")


def _write_bam(path, *, ref_name="REF", ref_len=10, n_primary=5, include_secondary=True):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": ref_name, "LN": ref_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i in range(n_primary):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * ref_len
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 60
            a.cigarstring = f"{ref_len}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * ref_len)
            bam.write(a)

        if include_secondary:
            sec = pysam.AlignedSegment()
            sec.query_name = "secondary"
            sec.query_sequence = "A" * ref_len
            sec.flag = 256  # secondary
            sec.reference_id = 0
            sec.reference_start = 0
            sec.mapping_quality = 60
            sec.cigarstring = f"{ref_len}M"
            sec.query_qualities = pysam.qualitystring_to_array("I" * ref_len)
            bam.write(sec)

    pysam.index(str(path))


def test_coverage_depth_counts_primary_only(tmp_path):
    from aavolve.coverage_depth import main

    bam = tmp_path / "reads.bam"
    _write_bam(bam, ref_len=10, n_primary=5, include_secondary=True)

    out = tmp_path / "depth.tsv.gz"
    rc = main(["--bam", str(bam), "--output", str(out)])
    assert rc == 0

    with gzip.open(out, "rt") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    assert lines[0] == "ref\tpos\tdepth"
    assert len(lines) == 1 + 10
    for row in lines[1:]:
        ref, pos, depth = row.split("\t")
        assert ref == "REF"
        assert int(pos) >= 1
        assert int(depth) == 5

