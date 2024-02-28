import tempfile
import sys

import pytest
import pysam

from scripts.extract_features_from_sam import (
    parse_args, get_reference_names, count_total_reads, check_index_exists,
    get_variants, identify_aa_change, get_query_base, write_header, write_variant,
    get_all_variants
)

from scripts.utils import Substitution, Deletion, Insertion

@pytest.fixture
def samfile_non_existent():
    return "non_existent.bam"

@pytest.fixture
def samfile_aav2():
    return "tests/data/aln/aav2_N496D.bam"

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

# ## Test parse_args
def test_parse_args_1():
    """
    Test parse_args with no arguments
    """
    with pytest.raises(SystemExit):
        parse_args([])

def test_parse_args_2():
    """
    Test parse_args with arguments
    """
    args = parse_args(["-i", "test.sam", "-r", "test.fasta", "-o", "test.tsv"])
    assert args.i == "test.sam"
    assert args.r == "test.fasta"
    assert args.o == "test.tsv"
    assert args.must_start_before == 0
    assert args.must_end_after == -1
    assert args.smaller_output == False

def test_parse_args_3():
    """
    Test parse_args with arguments
    """
    args = parse_args(["-i", "test.sam", "-r", "test.fasta", "-o", "test.tsv",
                       "--must-start-before", "1", "--must-end-after", "2",
                       "--smaller-output"])
    assert args.i == "test.sam"
    assert args.r == "test.fasta"
    assert args.o == "test.tsv"
    assert args.must_start_before == 1
    assert args.must_end_after == 2
    assert args.smaller_output == True

## Test check_index_exists
@pytest.mark.parametrize("samfile", ['samfile_aav2', 'samfile_pb', 'samfile_pb_sup_sec'], indirect=True)
def test_check_index_exists_1(samfile):
    """
    Test check_index_exists with small sam file
    """
    
    check_index_exists(samfile) 

def test_check_index_exists_2(samfile_non_existent):
    """
    Test check_index_exists with non-existent file
    """
    with pytest.raises(OSError):
        check_index_exists(samfile_non_existent)


### Test get_reference_names
def test_get_reference_names_1(samfile_aav2):
    """
    Test get_reference_names with small sam file
    """
    ref_names = get_reference_names(samfile_aav2)
    assert len(ref_names) == 1
    assert ref_names["AAV2_N496D"] == "AAV2"

@pytest.mark.parametrize("samfile", ['samfile_pb', 'samfile_pb_sup_sec'], indirect=True)
def test_get_reference_names_2(samfile):
    """
    Test get_reference_names with small sam file
    """
    ref_names = get_reference_names(samfile)
    assert len(ref_names) == 5
    assert all([x=="AAV2" for x in ref_names.values()])
    assert 'm54079_200305_102420/4194371/ccs' in ref_names.keys()
    assert 'm54079_200305_102420/4194437/ccs' in ref_names.keys()
    assert 'm54079_200305_102420/4194438/ccs' in ref_names.keys()
    assert 'm54079_200305_102420/4194440/ccs' in ref_names.keys()
    assert 'm54079_200305_102420/4194500/ccs' in ref_names.keys()


## Test count_total_reads
@pytest.mark.parametrize("alignmentfile,count", 
                         [('samfile_aav2', 1), 
                          ('samfile_pb', 5),
                          ('samfile_pb_sup_sec', 7)], 
                          indirect=['alignmentfile'])
def test_count_total_reads_1(alignmentfile, count):
    """
    Test count_total_reads with small sam file
    """
    assert count_total_reads(alignmentfile) == count

## Test get_variants
def test_get_variants_1(samfile_aav2, aav2_ref):
    """
    Test get_variants with small sam file
    """
    for read_id, var in get_variants(samfile_aav2, aav2_ref, start=0, end=-1, aa_isolation=False):
        assert read_id == "AAV2_N496D"
        assert len(var) == 1
        assert str(var[0]) == "A1486G"

@pytest.mark.parametrize("samfile,aav2_ref,start,end,count", [('samfile_pb','aav2_ref',0,-1,1), ('samfile_pb_sup_sec','aav2_ref',0,-1,1,)], indirect=['samfile', 'aav2_ref'])
def test_get_variants_2(samfile, aav2_ref, start, end, count):
    """
    Test get_variants with small sam file with many variants aav2_ref
    """
    variants = list(get_variants(samfile, aav2_ref, start=0, end=-1, aa_isolation=False))

    # check read names
    reads = [i[0] for i in variants]
    assert reads == ['m54079_200305_102420/4194371/ccs', 'm54079_200305_102420/4194437/ccs', 'm54079_200305_102420/4194438/ccs', 'm54079_200305_102420/4194440/ccs', 'm54079_200305_102420/4194500/ccs']

    # check number of variants
    vars = [len(i[1]) for i in variants]
    assert vars == [102, 75, 92, 93, 93]

    # first variant in each read is C403G
    assert all(str(i[1][0]) == 'C403G' for i in variants)


@pytest.mark.parametrize("start,end,samfile_toy,toy_ref", ((0,-1,"samfile_toy", "toy_ref"), (10, 90,"samfile_toy", "toy_ref")), indirect=['samfile_toy', 'toy_ref'])
def test_get_variants_3(start, end, samfile_toy, toy_ref):
    """
    Test get_variants with small sam file
    """
    variants = list(get_variants(samfile_toy, toy_ref, start=start, end=end, aa_isolation=False))

    if start == 0:
        assert len(variants) == 12
    else:
        assert len(variants) == 16

    # check each read

    # read 1 - matches reference
    assert variants[0] == ('read1_matches_ref', [])

    # read 2 - substitution
    assert variants[1][0] == 'read2_sub_G6A'
    assert [str(i) for i in variants[1][1]] == ['G6A']
    assert variants[1][1][0].changes_aa == True

    # read 3 - substitution and insertion before reference start (which isn't identified)
    assert variants[2][0] == 'read3_sub_A6G_insbeforeStart'
    assert [str(i) for i in variants[2][1]] == ['G6A']
    assert variants[2][1][0].changes_aa == True

    # read 4 - substitution and insertion after reference end 
    assert variants[3][0] == 'read4_sub_A6G_insafterend'
    assert [str(i) for i in variants[3][1]] == ['G6A', '100_101insTTTTT']
    assert variants[3][1][0].changes_aa == True

    # read 5 - insertion
    assert variants[4][0] == 'read5_7_8insttt'
    assert [str(i) for i in variants[4][1]] == ['7_8insTTT']

    # read 6 - deletion
    assert variants[5][0] == 'read6_T8_T10del'
    assert [str(i) for i in variants[5][1]] == ['T8_T10del']

    # read 7 - deletion and insertion
    assert variants[6][0] == 'read7_10_11_insaaaa_T18_T22_del'
    assert [str(i) for i in variants[6][1]] == ['10_11insAAAA', 'T18_T22del']

    # read 8 - deletion, insertion and substitution
    assert variants[7][0] == 'read8_10_11_insaaaa_T18_T22_del_subC29T'
    assert [str(i) for i in variants[7][1]] == ['10_11insAAAA', 'T18_T22del', 'C29T']
    assert variants[7][1][2].changes_aa == True

    # read 9 - sub first base
    assert variants[8][0] == 'read9_sub_C1G'
    assert [str(i) for i in variants[8][1]] == ['C1G']
    assert variants[8][1][0].changes_aa == True

    # read 10 - deletion of one base
    assert variants[9][0] == 'read10_T10del'
    assert [str(i) for i in variants[9][1]] == ['T10del']

    # read 11 - long insertion
    assert variants[10][0] == 'read11_50_51_insgtgtgtagagatttagagcgcgcttctcggatatatagcgcgcgcgttagag'
    assert [str(i) for i in variants[10][1]] == ['50_51insGTGTGTAGAGATTTAGAGCGCGCTTCTCGGATATATAGCGCGCGCGTTAGAG']

    # read 12 - long deletion
    assert variants[11][0] == 'read12_A28_A57del'
    assert [str(i) for i in variants[11][1]] == ['A28_A57del']

    # reads with deletions at start and end
    # these aren't in the CIGAR string, so are not identified
    # if start == 0 and end ==-1, these reads are skipped because the 
    # alignment doesn't cover the full reference
    # otherwise they're included
    if start != 0:
        # read 13 - deletion at end
        assert variants[12][0] == 'read13_sub_C1G_A98_A100del'
        assert [str(i) for i in variants[12][1]] == ['C1G']
        assert variants[12][1][0].changes_aa == True

        # read 14 - deletion at start
        assert variants[13][0] == 'read14_C1_T5del'
        assert [str(i) for i in variants[13][1]] == []

        # read 15 - deletion at start and end
        assert variants[14][0] == 'read15_C1_T5del_A98_A100del'
        assert [str(i) for i in variants[14][1]] == []

        # read 16 - deletion before start, insertion at start, deletion at end
        # deletions before start of reference are ignored
        assert variants[15][0] == 'read16_1_insaaa_C1_T5del_A98_A100del'
        assert [str(i) for i in variants[15][1]] == []

def check_var_type(read, var, aa_isolation):

    if var.var_type != 'sub':
        return True
    if 'nonsyn' in read:
        return True
    if not aa_isolation and 'synAndNonsyn' in read:
        return True
    if aa_isolation and 'synAndNonsyn' in read and str(var) == 'A62C':
        return True
    if aa_isolation and 'synAndNonsyn' in read and str(var) == 'G63A':
        return False
    if 'syn' in read:
        return False
    
    return None
    
        
@pytest.mark.parametrize('samfile_aav2_subs,aav2_ref,start,end,aa_isolation', 
                         (("samfile_aav2_subs", "aav2_ref", 0, -1, True),
                          ("samfile_aav2_subs", "aav2_ref", 0, -1, False),
                          ("samfile_aav2_subs", "aav2_ref", 5, -1, True),
                          ("samfile_aav2_subs", "aav2_ref", 5, -1, False),
                          ), 
                         indirect=['samfile_aav2_subs', 'aav2_ref'])
def test_get_variants_subs(samfile_aav2_subs, aav2_ref, start, end, aa_isolation):

    variants = list(get_variants(samfile_aav2_subs, aav2_ref, start=start, end=end, aa_isolation=aa_isolation))

    # check number of variants
    if start == 0:
        assert len(variants) == 11
    else:
        assert len(variants) == 12

    for read, vars in variants:
        for var in vars:
            try:
                assert var.changes_aa == check_var_type(read, var, aa_isolation)
            except AssertionError:
                import pdb; pdb.set_trace()
            

## Test identify_aa_change

@pytest.fixture
def mock_read(samfile_aav2):

    align = pysam.AlignmentFile(samfile_aav2, "rb")
    read = pysam.AlignedSegment(header=align.header)
    read.reference_name = "AAV2"
    return read

@pytest.mark.parametrize("aav2_ref,mock_read,query_seq,qpos,rpos,offset,aa_isolation,expected", 
                         (("aav2_ref", "mock_read",'ATG', 0, 0, 0, False, False), # no change
                          ("aav2_ref", "mock_read",'ATG', 1, 0, 0, False, False), # no change
                          ("aav2_ref", "mock_read",'ATGGCG', 5, 5, 0, False, False), # synonymous change
                          ("aav2_ref", "mock_read",'ATGCCG', 5, 5, 0, False, True), # different aa
                          ("aav2_ref", "mock_read",'ATGCG', 4, 5, -1, False, False), # synonymous change after deletion
                          ("aav2_ref", "mock_read",'ATGaGCG', 6, 5, 1, False, False), # synonymous change after insertion
                          ("aav2_ref", "mock_read",'ATGGAC', 5, 5, 0, True, False), # synonymous and non-synonymous changes - synonymous is specified, check in isolation
                          ("aav2_ref", "mock_read",'ATGGAC', 5, 5, 0, False, True), # synonymous and non-synonymous changes - synonymous is specified, don't check in isolation
                          ("aav2_ref", "mock_read",'ATGGAC', 4, 4, 0, True, True), # synonymous and non-synonymous changes - non-synonymous is specified, check in isolation
                          ("aav2_ref", "mock_read",'ATGGAC', 4, 4, 0, False, True), # synonymous and non-synonymous changes - non-synonymous is specified, don't check in isolation
                          ), 
                         indirect=['aav2_ref', 'mock_read'])
def test_identify_aa_change_1(aav2_ref, mock_read, query_seq, qpos, rpos, offset, aa_isolation, expected):
    """
    Test identify_aa_change with matching codons
    """
    
    mock_read.query_sequence = query_seq

    assert identify_aa_change(mock_read, aav2_ref, qpos, rpos, offset, aa_isolation) is expected


## Test get_query_base
@pytest.mark.parametrize("mock_read,query_seq,qpos,expected", 
                          (
                                ("mock_read", 'ATG', 0, 'A'),
                                ("mock_read", 'ATG', 1, 'T'),
                                ("mock_read", 'ATG', 2, 'G'),
                          ), 
                          indirect=['mock_read'])
def test_get_query_base_1(mock_read, query_seq, qpos,  expected):
    """
    Test get_query_base
    """
    
    mock_read.query_sequence = query_seq

    assert get_query_base(mock_read, qpos) == expected

## Test write_header

@pytest.mark.parametrize("smaller", (True, False))
def test_write_header(smaller):
    """
    Test write_header 
    """
    with tempfile.NamedTemporaryFile(mode='w+t') as filehandle:
        #import pdb; pdb.set_trace()
        write_header(filehandle.file, smaller=smaller)

        # read header
        filehandle.file.seek(0)
        header = filehandle.file.readlines()
    
    assert len(header) == 1
    
    header = header[0].strip().split("\t")
    if smaller:
        expected_header = ["query_name", "pos", "ref_bases", "query_bases", "aa_change"]
    else:
        expected_header = ["reference_name", "pos", "query_name", "var", "ref_bases", "query_bases", "aa_change"]
    assert header == expected_header


## Test write_variant
@pytest.mark.parametrize("smaller", (True, False))
def test_write_variant_sub(smaller):
    """
    Test write_variant 
    """

    read_id = 'test_read'
    ref_names = {'test_read': 'AAV2'}

    ins = Insertion(9, 10, "C")
    ins.add_another_base("G")
    ins.add_another_base("T")

    del1 = Deletion(5, 7, "A")
    del1.add_another_base("C")
    del1.add_another_base("G")

    del2 = Deletion(10, 11, "C")

    vars = (
        Substitution(rpos = 5, rseq = 'A', qseq = 'C', changes_aa=False),
        ins, del1, del2
    )
    
    with tempfile.NamedTemporaryFile(mode='w+t') as filehandle:
        #import pdb; pdb.set_trace()
        write_variant(read_id, vars, ref_names, filehandle.file, smaller=smaller)

        # read header
        filehandle.file.seek(0)
        lines = filehandle.file.readlines()
    
    assert len(lines) == 4
    
    if smaller:
        expected_lines = (["test_read", "5", "A", "C", "False"],
                          ["test_read", "10", ".", "CGT", "True"],
                          ["test_read", "5_8", "ACG", ".", "True"],
                          ["test_read", "10_11", "C", ".", "True"],
        )
    else:
        expected_lines = (["AAV2", "5", "test_read", "A6C", "A", "C", "False"],
                          ["AAV2", "10", "test_read", "10_11insCGT", ".", "CGT", "True"],
                          ["AAV2", "5_8", "test_read", "A6_G8del", "ACG", ".", "True"],
                          ["AAV2", "10_11", "test_read", "C11del", "C", ".", "True"],
        )
    expected_lines = ['\t'.join(i) + '\n' for i in expected_lines]

    assert expected_lines == lines


@pytest.mark.parametrize("samfile,reffile,resultfile, start_before, end_after, aa_isolation",
                         (("samfile_aav2", "aav2_ref_file", "resultfile_aav2", 0, -1, True),
                          ("samfile_aav2", "aav2_ref_file", "resultfile_aav2", 10, 2000, True),
                           ("samfile_aav2", "aav2_ref_file", "resultfile_aav2", 10, 4000, True),
                          ),
                         indirect=['samfile', 'reffile', 'resultfile'])
@pytest.mark.parametrize("smaller_output", (True, False))
def test_get_all_variants(samfile, reffile, resultfile, start_before, end_after, aa_isolation, smaller_output):


    with tempfile.NamedTemporaryFile(mode='w+t') as outfile:
        
        get_all_variants(samfile, reffile, outfile.name, start_before, end_after, smaller_output, aa_isolation)

        # read results
        outfile.seek(0)
        results = outfile.file.readlines()


    # read expected results
    with open(resultfile, 'rt') as handle:
        expected_results = handle.readlines()
        # for smaller output, rearrange and drop columns
        new_expected_results = []
        if smaller_output:
            for line in expected_results:
                line = line.strip().split('\t')
                line = [line[2], line[1], line[4], line[5], line[6]]
                line = '\t'.join(line) + '\n'
                new_expected_results.append(line)
            expected_results = new_expected_results
    if end_after > 2208:
        expected_results = expected_results[0:1] 

    assert results == expected_results

        


     
