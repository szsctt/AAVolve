import tempfile
import pytest
from scripts.utils import use_open, get_repeats_from_r2c2_name, seq_generator


class TestUseOpen:
    
    def test_read_file(self, fasta_file, fasta_lines):

        with use_open(fasta_file.name, 'rt') as handle:
            result = handle.readlines()
        
        assert result == fasta_lines
        fasta_file.close()
    
    def test_read_file_gz(self, fasta_file_gz, fasta_lines):
        with use_open(fasta_file_gz.name, 'rt') as handle:
            result = handle.readlines()
        
        assert result == fasta_lines
        fasta_file_gz.close()


class TestGetRepeatsFromR2C2:
    
    def test_get_repeats(self, fasta_contents):
        
        expected_repeats = [1, 1, 1, 2, 3, 5]
        repeats = [get_repeats_from_r2c2_name(name) for name in fasta_contents]
        assert repeats == expected_repeats

    def test_get_repeats_invalid(self):
        with pytest.raises(ValueError):
            get_repeats_from_r2c2_name('invalid_name')


class TestSeqGenerator:

    @pytest.mark.parametrize('fasta_file', ['fasta_file', 'fasta_file_gz'], indirect=True)
    def test_seq_generator(self, fasta_file, fasta_contents):
        n_seqs = 0
        with open(fasta_file.name, 'r') as handle:
            for name, seq in seq_generator(handle):
                name = name[1:] # seq generator doesn't remove the '>'
                assert name in fasta_contents
                assert seq == fasta_contents[name]
                n_seqs +=1 
        
        # check number of sequences is correct
        assert n_seqs == len(fasta_contents)
        fasta_file.close()

    def test_seq_generator_empty(self):
        # check empty file works
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            with open(temp.name, 'r') as handle:
                assert list(seq_generator(handle)) == []


class TestSubstituion:
    pass


class TestInsertion:
    pass

class TestDeletion:
    pass