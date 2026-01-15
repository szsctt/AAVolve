import gzip
import tempfile

import pytest

from aavolve.normalize_splint_fasta import main, normalize_splint_fasta


class TestNormalizeSplintFasta:
    def test_success_single_sequence(self):
        with tempfile.NamedTemporaryFile(mode="w+t") as in_fa, tempfile.NamedTemporaryFile(
            mode="w+t"
        ) as out_fa:
            in_fa.write(">original_name\nACGTACGT\n")
            in_fa.seek(0)

            normalize_splint_fasta(in_fa.name, out_fa.name)

            out_fa.seek(0)
            lines = out_fa.read().splitlines()

        assert lines[0] == ">splint"
        assert "".join(lines[1:]) == "ACGTACGT"

    def test_error_zero_sequences(self):
        with tempfile.NamedTemporaryFile(mode="w+t") as in_fa, tempfile.NamedTemporaryFile(
            mode="w+t"
        ) as out_fa:
            in_fa.write("")
            in_fa.seek(0)

            with pytest.raises(ValueError, match=r"exactly 1 sequence.*found 0"):
                normalize_splint_fasta(in_fa.name, out_fa.name)

    def test_error_multiple_sequences(self):
        with tempfile.NamedTemporaryFile(mode="w+t") as in_fa, tempfile.NamedTemporaryFile(
            mode="w+t"
        ) as out_fa:
            in_fa.write(">a\nAAAA\n>b\nCCCC\n")
            in_fa.seek(0)

            with pytest.raises(ValueError, match=r"exactly 1 sequence.*found 2"):
                normalize_splint_fasta(in_fa.name, out_fa.name)

    def test_gzipped_input(self):
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".gz") as in_gz, tempfile.NamedTemporaryFile(
            mode="w+t"
        ) as out_fa:
            with gzip.open(in_gz.name, "wt") as handle:
                handle.write(">x\nACGT\n")

            normalize_splint_fasta(in_gz.name, out_fa.name)

            out_fa.seek(0)
            lines = out_fa.read().splitlines()

        assert lines[0] == ">splint"
        assert "".join(lines[1:]) == "ACGT"

    def test_line_wrapping_60_chars(self):
        seq = "A" * 121
        with tempfile.NamedTemporaryFile(mode="w+t") as in_fa, tempfile.NamedTemporaryFile(
            mode="w+t"
        ) as out_fa:
            in_fa.write(">x\n" + seq + "\n")
            in_fa.seek(0)

            normalize_splint_fasta(in_fa.name, out_fa.name)

            out_fa.seek(0)
            lines = out_fa.read().splitlines()

        assert lines[0] == ">splint"
        assert len(lines[1]) == 60
        assert len(lines[2]) == 60
        assert len(lines[3]) == 1
        assert "".join(lines[1:]) == seq


class TestMain:
    def test_main_success(self):
        with tempfile.NamedTemporaryFile(mode="w+t") as in_fa, tempfile.NamedTemporaryFile(
            mode="w+t"
        ) as out_fa:
            in_fa.write(">x\nACGT\n")
            in_fa.seek(0)

            main(["--input", in_fa.name, "--output", out_fa.name])

            out_fa.seek(0)
            assert out_fa.read().splitlines()[0] == ">splint"
