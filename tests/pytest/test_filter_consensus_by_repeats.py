import tempfile
import gzip
from sys import argv

import pytest

from scripts.filter_consensus_by_repeats import seq_generator, filter, main

@pytest.fixture
def fasta_contents():
    seqs = {
        '11f28eab-af0c-4dee-94dd-57b84d8bc619_28.5_4474_1_2350': 'AGATAGATAGATGA',
        '8f7f7644-8296-43a5-9189-bc8ab62fc5c7_27.6_3179_1_2341': 'AGTAGCGAGCTACGAGAGCTATCGGCA',
        'f558aad8-dd60-44a2-86a0-6f9cd157e12b_20.5_5559_1_2354': 'AGATAGGAGATAGGAGCGCGTAGAGC',
        '399f92fe-a8b3-42f0-8515-7ae75d9ba16a_25.8_2856_2_2336': 'AGATAGAGCTGAGATCGCGGCTGG',
        '290de9a8-fd77-4c26-a56a-29799bbca0c7_23.7_5687_3_2356': 'AGTAGAGATCGCGTAGAGGTCGA',
        '39369093-a13a-43de-97f7-85cec3132850_24.5_2688_5_2336': 'AGATAGAGCTGAGATCGCGGCTGG',
    }
    return seqs

@pytest.fixture
def fasta_file(fasta_contents):
    temp = tempfile.NamedTemporaryFile(mode='w+t')
    for name, seq in fasta_contents.items():
        temp.write(f'>{name}\n{seq}\n')
    temp.seek(0)
    return temp

@pytest.fixture
def fasta_file_gz(fasta_contents):
    temp = tempfile.NamedTemporaryFile(mode='w+t', suffix='.gz')
    for name, seq in fasta_contents.items():
        temp.write(f'>{name}\n{seq}\n')
    temp.seek(0)
    return temp

class TestSeqGenerator:

    @pytest.mark.parametrize('fasta_file', [fasta_file, fasta_file_gz], indirect=True)
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

class TestFilter:

    @pytest.mark.parametrize('fasta_file', [fasta_file, fasta_file_gz], indirect=True)
    @pytest.mark.parametrize('min_repeats', [1, 2, 3, 4, 5])
    @pytest.mark.parametrize('out_gzipped', [True, False])
    def test_filter(self, fasta_file, min_repeats, out_gzipped):
        if out_gzipped:
            out_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.gz')
            out_open = gzip.open
        else:
            out_file = tempfile.NamedTemporaryFile(mode='w+t')
            out_open = open
        
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

    @pytest.mark.parametrize('fasta_file', [fasta_file, fasta_file_gz], indirect=True)
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
