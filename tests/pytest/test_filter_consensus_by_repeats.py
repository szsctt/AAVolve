import tempfile
import gzip
from sys import argv

import pytest

from scripts.filter_consensus_by_repeats import seq_generator, filter, main
from scripts.utils import use_open

class TestFilter:

    @pytest.mark.parametrize('fasta_file', ['fasta_file', 'fasta_file_gz'], indirect=True)
    @pytest.mark.parametrize('min_repeats', [1, 2, 3, 4, 5])
    @pytest.mark.parametrize('out_gzipped', [True, False])
    def test_filter(self, fasta_file, min_repeats, out_gzipped):
        # create output file with correct file extension
        if out_gzipped:
            out_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.gz')
        else:
            out_file = tempfile.NamedTemporaryFile(mode='w+t')
        
        # call function
        filter(fasta_file.name, out_file.name, min_repeats)

        # check output
        out_file.seek(0)
        with use_open(out_file.name, 'rt') as handle:
            for line in handle:
                if line.startswith('>'):
                    name = line.strip()
                    repeats = int(name.split('_')[-2])
                    assert repeats >= min_repeats
        
        fasta_file.close()
        out_file.close()

    def test_filter_empty(self):
        # check empty file works
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            with tempfile.NamedTemporaryFile(mode='w+t') as out_temp:
                filter(temp.name, out_temp.name, 1)
                out_temp.seek(0)
                assert out_temp.read() == ''

    def test_filter_invalid(self):
        # check invalid min_repeats raises error
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            with tempfile.NamedTemporaryFile(mode='w+t') as out_temp:
                with pytest.raises(ValueError):
                    filter(temp.name, out_temp.name, 0)

    @pytest.mark.parametrize('out_gzipped', [True,False])
    def test_filter_invalid_outfile(self, out_gzipped):
        if out_gzipped:
            out_file = 'test.gz'
        else:
            out_file = 'test'
        with pytest.raises(FileNotFoundError):
            filter('doesnotexist', out_file, 1)

class TestMain:

    @pytest.mark.parametrize('fasta_file', ['fasta_file', 'fasta_file_gz'], indirect=True)
    @pytest.mark.parametrize('min_repeats', [1, 2, 3, 4, 5])
    @pytest.mark.parametrize('out_gzipped', [True, False])
    def test_main(self, fasta_file, min_repeats, out_gzipped):
        
        if out_gzipped:
            out_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.gz')
            out_open = gzip.open
        else:
            out_file = tempfile.NamedTemporaryFile(mode='w+t')
            out_open = open
        
        main(['--input', fasta_file.name, '--output', out_file.name, '--min-repeats', '2'])
        
        filter(fasta_file.name, out_file.name, min_repeats)
        out_file.seek(0)
        with out_open(out_file.name, 'rt') as handle:
            for line in handle:
                if line.startswith('>'):
                    name = line.strip()
                    repeats = int(name.split('_')[-2])
                    assert repeats >= min_repeats
        
        fasta_file.close()
        out_file.close()

    def test_main_invalid_min_repeats(self):
        # check invalid min_repeats raises error
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            with tempfile.NamedTemporaryFile(mode='w+t') as out_temp:
                with pytest.raises(ValueError):
                    main(['--input', temp.name, '--output', out_temp.name, '--min-repeats', '0'])
