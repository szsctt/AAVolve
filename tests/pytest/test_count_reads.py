import tempfile
import gzip
import pytest

from scripts.count_reads import get_file_type, count_fasta, count_fastq, count_variant_tsv, count_pivoted_tsv, main

@pytest.fixture
def write_fa_no_suffix():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        f.write('>seq1\nACGT\n>seq2\nACGT\n')
        f.seek(0)
        yield f.name

@pytest.fixture
def write_fna():
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.fna') as f:
        f.write('>seq1\nACGT\n>seq2\nACGT\n')
        f.seek(0)
        yield f.name    

@pytest.fixture
def empty_fa():
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.fa') as f:
        yield f.name   

@pytest.fixture
def np_2389_reads():
    return "tests/data/reads/np-2389.fastq"

@pytest.fixture
def np_cc_aav2_reads():
    return "tests/data/reads/np-cc-aav2.fastq"

def gzip_file(filename):
    g = tempfile.NamedTemporaryFile(mode='w+t', suffix = ".gz")
    with open(filename, 'rt') as f:
        with gzip.open(g.name, 'wt') as h:
            h.write(f.read())
    return g

@pytest.fixture
def gzipped_reads(np_2389_reads):
    
    gz = gzip_file(np_2389_reads)
    yield gz.name
    gz.close()

@pytest.fixture
def empty_fq():
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.fq') as f:
        yield f.name  

@pytest.fixture
def gzipped_variants(resultfile_aav2):
    gz = gzip_file(resultfile_aav2)
    yield gz.name
    gz.close()

@pytest.fixture
def empty_variants():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        # write header
        f.write('reference_name\tpos\tquery_name\tvar\tref_base\tquery_bases\taa_change\n')
        f.seek(0)
        yield f.name

@pytest.fixture
def gzipped_pivoted_variants(aav3_pivoted_seq):
    gz = gzip_file(aav3_pivoted_seq)
    yield gz.name
    gz.close()

@pytest.fixture
def empty_pivoted():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        # write header
        f.write('read_id\tsub:40\tsub:50\n')
        f.seek(0)
        yield f.name

@pytest.fixture
def empty_pivoted_no_vars():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        # write header
        f.write('read_id\n')
        f.seek(0)
        yield f.name

@pytest.fixture
def gzipped_counts_file(counts_file):
    gz = gzip_file(counts_file)
    yield gz.name
    gz.close()

@pytest.fixture
def empty_file():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        yield f.name     

@pytest.fixture
def input_file(request):
    return request.getfixturevalue(request.param)

class TestGetFileType:

    @pytest.mark.parametrize("input_file,expected", (
            ('aav2_ref_file', "fasta"),
            ('aav2389_ref_file', "fasta"),
            ('write_fa_no_suffix', "fasta"),
            ('write_fna', "fasta"),
            ('empty_fa', 'fasta'),
            ('np_2389_reads', "fastq"),
            ('np_cc_aav2_reads', "fastq"),
            ('gzipped_reads', "fastq"),
            ('empty_fq', 'fastq'),
            ('resultfile_aav2', 'variant_tsv'),
            ('resultfile_aav2_shorter', 'variant_tsv'),
            ('gzipped_variants', 'variant_tsv'),
            ('empty_variants', 'variant_tsv'),
            ('aav3_pivoted_seq', 'pivoted_tsv'),
            ('gzipped_pivoted_variants', 'pivoted_tsv'),
            ('empty_pivoted', 'pivoted_tsv'),
            ('empty_pivoted_no_vars', 'pivoted_tsv'),
            ('counts_file', 'distinct_read_counts'),
            ('gzipped_counts_file', 'distinct_read_counts'),
            ('empty_counts', 'distinct_read_counts'),
            ('empty_file', 'unknown'),
            ('samfile_toy', 'unknown'),

    ), indirect=['input_file'])
    def test_get_file_type(self, input_file, expected):
        
        assert get_file_type(input_file) == expected

class TestCountFasta:

    @pytest.mark.parametrize("input_file,expected", (
            ('aav2_ref_file', 1),
            ('aav2389_ref_file', 4),
            ('write_fa_no_suffix', 2),
            ('write_fna', 2),
            ('empty_fa', 0),
            ('np_2389_reads', 0),
            ('np_cc_aav2_reads', 0),
            ('gzipped_reads', 0),
            ('empty_fq', 0),
            ('resultfile_aav2', 0),
            ('resultfile_aav2_shorter', 0),
            ('gzipped_variants', 0),
            ('empty_variants', 0),
            ('aav3_pivoted_seq', 0),
            ('gzipped_pivoted_variants', 0),
            ('empty_pivoted', 0),
            ('empty_pivoted_no_vars', 0),
            ('counts_file', 0),
            ('gzipped_counts_file', 0),
            ('empty_counts', 0),
            ('empty_file', 0),
            ('samfile_toy', 0),

    ), indirect=['input_file'])
    def test_count_fasta(self, input_file, expected):
        
        assert count_fasta(input_file) == expected

class TestCountFastq:

    @pytest.mark.parametrize("input_file,expected", (
            ('aav2_ref_file', 0),
            ('aav2389_ref_file', 0),
            ('write_fa_no_suffix', 0),
            ('write_fna', 0),
            ('empty_fa', 0),
            ('np_2389_reads', 1000),
            ('np_cc_aav2_reads', 250),
            ('gzipped_reads', 1000),
            ('empty_fq', 0),
            ('resultfile_aav2', 0),
            ('resultfile_aav2_shorter', 0),
            ('gzipped_variants', 0),
            ('empty_variants', 0),
            ('aav3_pivoted_seq', 0),
            ('gzipped_pivoted_variants', 0),
            ('empty_pivoted', 0),
            ('empty_pivoted_no_vars', 0),
            ('counts_file', 0),
            ('gzipped_counts_file', 0),
            ('empty_counts', 0),
            ('empty_file', 0),
            ('samfile_toy', 0),

    ), indirect=['input_file'])
    def test_count_fastq(self, input_file, expected):
        
        assert count_fastq(input_file) == expected

class TestCountVariantTsv:

    @pytest.mark.parametrize("input_file,expected", (
            ('resultfile_aav2', 1),
            ('resultfile_aav2_shorter', 1),
            ('gzipped_variants', 1),
            ('resultfile_aav23', 1),
            ('resultfile_toy', 11),
            ('empty_variants', 0),

    ), indirect=['input_file'])
    def test_count_variant_tsv(self, input_file, expected):
        
        assert count_variant_tsv(input_file) == expected

class TestCountPivotedTsv:

    @pytest.mark.parametrize("input_file,expected", (
            ('aav3_pivoted_seq', 1),
            ('gzipped_pivoted_variants', 1),
            ('toy_pivoted_seq', 11),
            ('toy_pivoted_parents', 11),
            ('empty_pivoted', 0),
            ('empty_pivoted_no_vars', 0),
            ('counts_file', 2),
            ('gzipped_counts_file', 2),
            ('empty_counts', 0),
    ), indirect=['input_file'])
    def test_count_pivoted_tsv(self, input_file, expected):
        
        assert count_pivoted_tsv(input_file) == expected

class TestMain:

    @pytest.mark.parametrize("input_file,expected_file_type,expected_count", (
            ('aav2_ref_file', "fasta", 1),
            ('aav2389_ref_file', "fasta", 4),
            ('write_fa_no_suffix', "fasta", 2),
            ('write_fna', "fasta", 2),
            ('empty_fa', 'fasta', 0),
            ('np_2389_reads', "fastq", 1000),
            ('np_cc_aav2_reads', "fastq", 250),
            ('gzipped_reads', "fastq", 1000),
            ('empty_fq', 'fastq', 0),
            ('resultfile_aav2', 'variant_tsv', 1),
            ('resultfile_aav2_shorter', 'variant_tsv', 1),
            ('gzipped_variants', 'variant_tsv', 1),
            ('empty_variants', 'variant_tsv', 0),
            ('aav3_pivoted_seq', 'pivoted_tsv', 1),
            ('gzipped_pivoted_variants', 'pivoted_tsv', 1),
            ('empty_pivoted', 'pivoted_tsv', 0),
            ('empty_pivoted_no_vars', 'pivoted_tsv', 0),
            ('counts_file', 'distinct_read_counts', 2),
            ('gzipped_counts_file', 'distinct_read_counts', 2),
            ('empty_counts', 'distinct_read_counts', 0),
            ('empty_file', 'unknown', 0),
            ('samfile_toy', 'unknown', 0),
    ), indirect=['input_file'])
    def test_main(self, input_file, expected_file_type, expected_count, monkeypatch):
        
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            # set sys.argv
            monkeypatch.setattr('sys.argv', ['count_reads.py', '--output', f.name, '--files', input_file])

            # run main
            if expected_file_type == "unknown":
                with pytest.raises(ValueError):
                    main()
                return
            main()

            # check output
            f.seek(0)
            assert f.read() == f"{input_file}\t{expected_file_type}\t{expected_count}\n"

    def test_main_multiple(self, write_fna, np_2389_reads, resultfile_aav2, aav3_pivoted_seq, counts_file, monkeypatch):
        
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            # set sys.argv
            monkeypatch.setattr('sys.argv', ['count_reads.py', '--output', f.name, '--files', write_fna, np_2389_reads, resultfile_aav2, aav3_pivoted_seq, counts_file])

            # run main
            main()

            # check output
            f.seek(0)
            assert f.read() == f"{write_fna}\tfasta\t2\n{np_2389_reads}\tfastq\t1000\n{resultfile_aav2}\tvariant_tsv\t1\n{aav3_pivoted_seq}\tpivoted_tsv\t1\n{counts_file}\tdistinct_read_counts\t2\n"


