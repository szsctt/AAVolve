import tempfile
import pytest

from scripts.num_in_fa import count_records, write_output


class TestCountRecords:

    @pytest.mark.parametrize('fasta_file', ['fasta_file', 'fasta_file_gz'], indirect=True)
    def test_count_records(self, fasta_file):
        expected_count = 6
        count = count_records(fasta_file.name)
        assert count == expected_count
        fasta_file.close()

    def test_count_records_empty(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            count = count_records(temp.name)
            assert count == 0

class TestWriteOutput:

    @pytest.mark.parametrize('n', [6, 0])
    def test_write_output(self, n):
        expected_result = f'{n}\n'

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            write_output(temp.name, n)
            temp.seek(0)
            result = temp.read()
        
        assert result == expected_result