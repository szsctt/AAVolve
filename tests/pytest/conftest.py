import tempfile
import gzip
import pytest
import pysam
import pandas as pd
from Bio import SeqIO

from aavolve.utils import Substitution, Insertion, Deletion
from aavolve.extract_features_from_sam import main as extract_features_from_sam
from aavolve.pivot_variants_to_wide import main as pivot_variants_to_wide


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
def aav2n496d_ref_file():
    return "tests/data/references/wtAAV2_N496D.fa"

@pytest.fixture
def aav2389_ref_file():
    return "tests/data/references/wt2n496d389dna.fasta"

@pytest.fixture
def aav23_ref_file():
    return "tests/data/references/AAV2_AAV3.fa"

@pytest.fixture
def toy_ref_file():
    return "tests/data/references/toy_reference.fa"

@pytest.fixture
def toy_reads_file():
    return "tests/data/reads/toy_reads.fa"

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


#### samfile fixtures ####
@pytest.fixture
def samfile_non_existent():
    return "non_existent.bam"

@pytest.fixture
def samfile_aav2_wt():
    return "tests/data/aln/aav2.bam"

@pytest.fixture
def samfile_aav2():
    return "tests/data/aln/aav2_N496D.bam"

@pytest.fixture
def samfile_aav23():
    return "tests/data/aln/aav23.bam"

@pytest.fixture
def samfile_aav2389():
    return "tests/data/aln/aav2389.bam"

@pytest.fixture
def alignmentfile_aav2(samfile_aav2):
    return pysam.AlignmentFile(samfile_aav2, "rb")

@pytest.fixture
def samfile_pb():
    return "tests/data/aln/pb-shuf.bam"

@pytest.fixture
def samfile_pb_sup_sec():
    return "tests/data/aln/pb-shuf_sup_sec.bam"

@pytest.fixture
def samfile_toy():
    return "tests/data/aln/toy.bam"

@pytest.fixture
def samfile_aav2_subs():
    return "tests/data/aln/aav2_subs.bam"

@pytest.fixture
def samfile(request):
    return request.getfixturevalue(request.param)  

@pytest.fixture
def alignmentfile(request):

    filename = request.getfixturevalue(request.param)
    return pysam.AlignmentFile(filename, "rb")

@pytest.fixture
def reffile(request):
    return request.getfixturevalue(request.param)


#### variant fixtures ####

def get_variants(sam, reference, *kwargs):

    f = tempfile.NamedTemporaryFile(mode='w+t')
    f2 = tempfile.NamedTemporaryFile(mode='w+t')
    extract_features_from_sam(['-i', sam, '-r', reference, '-o', f.name, '-O', f2.name, *kwargs])
    f.seek(0), f2.seek(0)
    return f, f2


@pytest.fixture()
def resultfile_aav2(aav2_ref_file, samfile_aav2):

    f, f2 = get_variants(samfile_aav2, aav2_ref_file)
    yield f.name, f2.name
    f.close(), f2.close()
    
@pytest.fixture
def resultfile_aav2_shorter(aav2_ref_file, samfile_aav2):

    f, f2 = get_variants(samfile_aav2, aav2_ref_file, '-S')
    yield f.name, f2.name
    f.close(), f2.close()

@pytest.fixture
def resultfile_aav23(aav2_ref_file, samfile_aav23):

    f, f2 = get_variants(samfile_aav23, aav2_ref_file)
    yield f.name, f2.name
    f.close(), f2.close()

@pytest.fixture
def resultfile_toy(toy_ref_file, samfile_toy):

    f, f2 = get_variants(samfile_toy, toy_ref_file)
    yield f.name, f2.name
    f.close(), f2.close()

@pytest.fixture
def resultfile_aav2389(aav2_ref_file, samfile_aav2389):
    
    f , f2 = get_variants(samfile_aav2389, aav2_ref_file)
    yield f.name, f2.name
    f.close(), f2.close()

@pytest.fixture
def resultfile_aav2389_some():
    return "tests/data/variants/test_variants.tsv", "tests/data/variants/test_variants_read_ids.tsv"

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
def resultfile_aav2389_some2_variants():
    return {
                '40:sub': {'AAV3b': Substitution(40, 'C', 'A', True),
                           'AAV9': Substitution(40, 'C', 'A', True),
                           'AAV8': Substitution(40, 'C', 'A', True),
                           },
                '44:sub': {'AAV3b': Substitution(44, 'C', 'T', False),
                           'AAV9': Substitution(44, 'C', 'T', False),
                           },
                '45:sub': {'AAV9': Substitution(45, 'T', 'A', False)
                           },
                '46:sub': {'AAV9': Substitution(46, 'C', 'G', False)
                           },
                '50:sub': {'AAV8': Substitution(50, 'A', 'G', False),
                           },
                '53:sub': {'AAV3b': Substitution(53, 'A', 'C', False),
                           'AAV8': Substitution(53, 'A', 'C', False),
                           },
        }

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

    with tempfile.NamedTemporaryFile(mode='w+t') as temp_vars, tempfile.NamedTemporaryFile(mode='w+t') as temp_rids:
        # write variants
        temp_vars.write(variants[0].header(shorter=shorter))
        for var in variants:
            reference_name = None if shorter else 'reference'
            temp_vars.write(var.print_line('query', reference_name))
        temp_vars.seek(0)
        # write read ids
        temp_rids.write("query\n")
        temp_rids.seek(0)

        yield shorter, variants, temp_vars, temp_rids

@pytest.fixture
def write_variants_repeated(request):

    shorter, n_repeats, variants_name = request.param
    variants = request.getfixturevalue(variants_name)
    assert len(n_repeats) == len(variants) # n_repeats for each variant

    with tempfile.NamedTemporaryFile(mode='w+t') as temp_vars, tempfile.NamedTemporaryFile(mode='w+t') as temp_rids:
        
        # write variants and read ids
        written_queries = set()
        temp_vars.write(variants[0].header(shorter=shorter))
        for i, var in enumerate(variants):
            reference_name = None if shorter else 'reference'
            for j in range(n_repeats[i]):
                temp_vars.write(var.print_line(f'query{j}', reference_name))
                if f'query{j}' not in written_queries:
                    temp_rids.write(f"query{j}\n")
                    written_queries.add(f'query{j}')
        temp_vars.seek(0)
        temp_rids.seek(0)

        yield shorter, n_repeats, variants, temp_vars, temp_rids

#### variant_pivoted fixtures ####

def pivot(variants, read_ids, parents):

    f_seq, f_parents = tempfile.NamedTemporaryFile(mode='w+t'), tempfile.NamedTemporaryFile(mode='w+t')
    pivot_variants_to_wide(['-i', variants, '-r', read_ids, '-p', parents, '-o', f_parents.name, '-O', f_seq.name])
    f_seq.seek(0), f_parents.seek(0)
    return f_seq, f_parents

@pytest.fixture
def aav3_pivoted_seq(resultfile_aav23):

    f_seq, f_parents = pivot(resultfile_aav23[0], resultfile_aav23[1], resultfile_aav23[0])
    yield f_seq.name
    f_seq.close(), f_parents.close()


@pytest.fixture
def aav3_pivoted_parents(resultfile_aav23):


    f_seq, f_parents = pivot(resultfile_aav23[0], resultfile_aav23[1], resultfile_aav23[0])
    yield f_parents.name
    f_seq.close(), f_parents.close()
    
@pytest.fixture
def toy_pivoted_seq(resultfile_toy):

    f_seq, f_parents = pivot(resultfile_toy[0], resultfile_toy[1], resultfile_toy[0])
    yield f_seq.name
    f_seq.close(), f_parents.close()

@pytest.fixture
def toy_pivoted_parents(resultfile_toy):
    f_seq, f_parents = pivot(resultfile_toy[0], resultfile_toy[1], resultfile_toy[0])
    yield f_parents.name
    f_seq.close(), f_parents.close()

#### counts ####
    
@pytest.fixture
def counts_file():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        f.write('count\tsequence\n1\tACGT\n2\tACGT\n')
        f.seek(0)
        yield f.name

@pytest.fixture
def empty_counts():
    with tempfile.NamedTemporaryFile(mode='w+t') as f:
        # write header
        f.write('count\tsequence\n')
        f.seek(0)
        yield f.name