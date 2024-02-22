import tempfile
import pytest

from scripts.first_last_variant import get_first_last_variant, write_output

class TestGetFirstLastVariant:

    @pytest.mark.parametrize("resultfile, expected_first, expected_last", [
        ("resultfile_aav2", 1485, 1485),
        ("resultfile_aav23", 8, 2206),
        ("resultfile_aav2389", 40, 2202),
    ], indirect=["resultfile"])
    def test_get_first_last_variant(self, resultfile, expected_first, expected_last):
        first, last = get_first_last_variant(resultfile)
        assert first == expected_first
        assert last == expected_last


    def test_get_first_last_variant_empty(self):

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            with pytest.raises(AssertionError):
                _ = get_first_last_variant(temp.name)

class TestWriteOutput:

    def test_write_output(self):

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            write_output(temp.name, 1, 2)
            temp.seek(0)
            assert temp.read() == "1\n2\n"


