import gzip

import pandas as pd

from aavolve.combine_non_parental_variants import combine_non_parental_variants
from aavolve.snakemake_helpers import build_non_parental_variant_group_maps
from aavolve.utils import read_variant_file, get_variant


def _write_gz(path: str, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def _variants_in_file(path: str) -> list[str]:
    return [str(get_variant(row)) for row in read_variant_file(path)]


def test_combine_non_parental_variants_includes_variants_from_other_samples(tmp_path):
    header = "reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n"

    sample_early = tmp_path / "early_high.tsv.gz"
    _write_gz(str(sample_early), header)

    sample_late = tmp_path / "late_high.tsv.gz"
    _write_gz(
        str(sample_late),
        header
        + "ref\t0\tnon_parental_0\tA1C\tA\tC\tFalse\n",
    )

    combined = tmp_path / "combined.tsv.gz"
    count = combine_non_parental_variants([str(sample_early), str(sample_late)], str(combined))

    assert count == 1
    assert _variants_in_file(str(combined)) == ["A1C"]


def test_combine_non_parental_variants_deduplicates_and_renumbers(tmp_path):
    header = "reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n"

    a = tmp_path / "a.tsv.gz"
    _write_gz(
        str(a),
        header
        + "ref\t0\tnon_parental_0\tA1C\tA\tC\tFalse\n"
        + "ref\t1\tnon_parental_1\tG2T\tG\tT\tFalse\n",
    )

    b = tmp_path / "b.tsv.gz"
    _write_gz(
        str(b),
        header
        + "ref\t0\tnon_parental_0\tA1C\tA\tC\tFalse\n"
        + "ref\t2\tnon_parental_1\tT3G\tT\tG\tFalse\n",
    )

    out = tmp_path / "out.tsv.gz"
    combine_non_parental_variants([str(a), str(b)], str(out))

    assert set(_variants_in_file(str(out))) == {"A1C", "G2T", "T3G"}

    rows = list(read_variant_file(str(out)))
    assert [row["query_name"] for row in rows] == ["non_parental_0", "non_parental_1", "non_parental_2"]


def test_combine_non_parental_variants_keeps_multiple_alleles_at_same_position(tmp_path):
    header = "reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n"

    a = tmp_path / "a.tsv.gz"
    _write_gz(str(a), header + "ref\t0\tnon_parental_0\tA1C\tA\tC\tFalse\n")

    b = tmp_path / "b.tsv.gz"
    _write_gz(str(b), header + "ref\t0\tnon_parental_0\tA1G\tA\tG\tFalse\n")

    out = tmp_path / "out.tsv.gz"
    combine_non_parental_variants([str(a), str(b)], str(out))

    assert set(_variants_in_file(str(out))) == {"A1C", "A1G"}


def test_build_non_parental_variant_group_maps(tmp_path):
    samples = pd.DataFrame(
        [
            {
                "sample_name": "S1",
                "include_non_parental": True,
            },
            {
                "sample_name": "S2",
                "include_non_parental": True,
            },
            {
                "sample_name": "S3",
                "include_non_parental": False,
            },
        ]
    )

    input_validation_samples = {"pair": ["S1", "S2", "S3"]}
    sample_to_pair_id, group_map = build_non_parental_variant_group_maps(samples, input_validation_samples)

    assert sample_to_pair_id["S1"] == sample_to_pair_id["S2"] == sample_to_pair_id["S3"] == "pair"
    assert list(group_map.values()) == [["S1", "S2"]]
