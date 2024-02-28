import tempfile
import gzip
import pytest
import pandas as pd
from Bio import SeqIO

from scripts.utils import Substitution, Insertion, Deletion

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

#### reference fasta fixtures ####

@pytest.fixture
def aav2_ref_file():
    return "tests/data/references/wtAAV2.fa"

@pytest.fixture
def aav2389_ref_file():
    return "tests/data/references/wt2n496d389dna.fasta"

@pytest.fixture
def aav2_ref(aav2_ref_file):
    return SeqIO.to_dict(SeqIO.parse(aav2_ref_file, "fasta"))

@pytest.fixture
def toy_ref():
    return SeqIO.to_dict(SeqIO.parse('tests/data//references/toy_reference.fa', "fasta"))

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
def resultfile_aav2389_some():
    return "tests/data/variants/test_variants.tsv"

@pytest.fixture
def resultfile_aav2389_some2():
    return "tests/data/variants/test_parents.tsv"

@pytest.fixture
def resultfile(request):    
    return request.getfixturevalue(request.param)


@pytest.fixture
def some_variants():
    vars = [
            Substitution(10, 'A', 'T', True),
            Substitution(13, 'A', 'G', False),
            Insertion(70, 82, 'C'),
            Insertion(75, 90, 'C'),
            Deletion(49, 61, 'C'),
            Deletion(90, 81, 'C')
        ]
    vars[2].add_another_base('G')
    vars[2].add_another_base('C')
    vars[5].add_another_base('G')
    vars[5].add_another_base('C')

    return vars

@pytest.fixture
def some_variant_frequencies(some_variants):
    freqs = dict()
    for i, var in enumerate(some_variants):
        freqs[var] = 1/(i+1)

    return freqs

@pytest.fixture
def write_vars(request):

    shorter, variants_name = request.param
    variants = request.getfixturevalue(variants_name)

    temp = tempfile.NamedTemporaryFile(mode='w+t')
    temp.write(variants[0].header(shorter=shorter))
    for var in variants:
        reference_name = None if shorter else 'reference'
        temp.write(var.print_line('query', reference_name))
    temp.seek(0)
    return shorter, variants, temp

@pytest.fixture
def write_variants_repeated(request):

    shorter, n_repeats, variants_name = request.param
    variants = request.getfixturevalue(variants_name)
    assert len(n_repeats) == len(variants) # n_repeats for each variant

    temp = tempfile.NamedTemporaryFile(mode='w+t')
    temp.write(variants[0].header(shorter=shorter))
    for i, var in enumerate(variants):
        reference_name = None if shorter else 'reference'
        for j in range(n_repeats[i]):
            temp.write(var.print_line(f'query{j}', reference_name))

    temp.seek(0)
    return shorter, n_repeats, variants, temp

#### variant_pivoted fixtures ####

def aav3_pivoted():
    'tests/data/variants/aav3_pivoted.tsv'

