import tempfile
import gzip
import pytest
import pandas as pd

#### sample csv fixtures ####

@pytest.fixture
def config():
    return {
        'sample_name': 'sample1',
        'parent_name': 'aavs',
        'reference_name': 'AAV2',
        'seq_tech': 'np',
        'read_file': 'tests/data/reads/np-aav2.fastq',
        'parent_file': 'tests/data/references/AAV2_AAV3.fa',
        'reference_file': 'tests/data/references/wtAAV2.fa'
    }

@pytest.fixture 
def sample_df(config):

    return pd.DataFrame([config])


#### nanopore RCA fasta fixtures ####

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
def fasta_lines(fasta_contents):
    contents = []
    for name, seq in fasta_contents.items():
        contents.append(f'>{name}\n')
        contents.append(f'{seq}\n')
    return contents

@pytest.fixture
def fasta_file(fasta_lines):
    temp = tempfile.NamedTemporaryFile(mode='w+t')
    for line in fasta_lines:
        temp.write(line)
    temp.seek(0)
    return temp

@pytest.fixture
def fasta_file_gz(fasta_lines):
    temp = tempfile.NamedTemporaryFile(mode='w+t', suffix='.gz')
    with gzip.open(temp.name, 'wt') as handle:
        for line in fasta_lines:
            handle.write(line)
    temp.seek(0)
    return temp


#### variant fixtures ####

@pytest.fixture
def resultfile_aav2():
    return "tests/data/variants/aav2N496D.tsv"

@pytest.fixture
def resultfile_aav2_shorter():
    return "tests/data/variants/aav2N496D_shorter.tsv"

@pytest.fixture
def resultfile_aav23():
    return "tests/data/variants/aav23.tsv.gz"

@pytest.fixture
def resultfile_aav2389():
    return "tests/data/variants/aav2389.tsv.gz"

@pytest.fixture
def resultfile(request):    
    return request.getfixturevalue(request.param)