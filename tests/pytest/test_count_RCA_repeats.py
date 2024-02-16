import tempfile
import pytest

from scripts.count_RCA_repeats import count_repeats, write_output

class TestCountRepeats:

    @pytest.mark.parametrize('fasta_file', ['fasta_file', 'fasta_file_gz'], indirect=True)
    def test_count_repeats(self, fasta_file):
        expected_repeats = {1: 3, 2: 1, 3: 1, 5: 1}
        repeats = count_repeats(fasta_file.name)
        assert repeats == expected_repeats
        fasta_file.close()

    def test_count_repeats_empty(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            repeats = count_repeats(temp.name)
            assert repeats == {}

    def test_count_repeats_non_RCA(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            temp.write('>seq1\nATGC\n')
            temp.seek(0)
            with pytest.raises(ValueError):
                repeats = count_repeats(temp.name)


class TestWriteOutput:
    def test_write_output(self):
        repeats = {1: 3, 2: 1, 3: 1, 5: 1}
        expected_result = 'Repeats\tCount\n1\t3\n2\t1\n3\t1\n5\t1\n'

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            write_output(temp.name, repeats)
            temp.seek(0)
            result = temp.read()
        
        assert result == expected_result