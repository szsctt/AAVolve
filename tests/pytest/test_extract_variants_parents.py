import pytest

from aavolve.utils import read_variant_file


def test_extract_variants_parents_all_parents_present(resultfile_aav2389_some2, resultfile_aav2389_some2_variants):
    """
    Ensure that the parents variants file contains variants for all expected parents.
    Uses the provided fixture `resultfile_aav2389_some2` which is a path to a
    variants TSV for parents, and `resultfile_aav2389_some2_variants` which
    describes the expected parent mapping.
    """

    variants_file = resultfile_aav2389_some2

    # read parents present in the variants file
    parents_in_file = set()
    for row in read_variant_file(variants_file):
        # long format has a 'query_name' field containing parent name
        if 'query_name' in row:
            parents_in_file.add(row['query_name'])

    # expected parents derived from the fixture mapping
    expected_parents = set()
    for var_map in resultfile_aav2389_some2_variants.values():
        expected_parents.update(var_map.keys())

    assert expected_parents, "Sanity check: expected parents should not be empty"
    assert parents_in_file == expected_parents
