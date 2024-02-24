import tempfile
import pytest
from scripts.utils import (
    use_open, get_repeats_from_r2c2_name, seq_generator, 
    read_variant_file, get_variant_type, get_variant, 
    get_variants_set,
    Substitution, Insertion, Deletion
    )

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

class TestReadVariantFile:

    @pytest.mark.parametrize('resultfile,expected_n_lines', 
                                [('resultfile_aav2',1), 
                                 ('resultfile_aav2_shorter',1), 
                                 ('resultfile_aav23', 416), 
                                 ('resultfile_aav2389', 1362)], 
                            indirect=['resultfile'])
    def test_read_variant_file(self, resultfile, expected_n_lines):
        n_lines = 0
        longer_cols = set(["reference_name", "pos", "query_name", "var", "ref_bases", "query_bases", "aa_change"])
        shorter_cols = set(["query_name", "pos", "ref_bases", "query_bases", "aa_change"])
        for line in read_variant_file(resultfile):
            assert line
            assert len(line) in [len(longer_cols), len(shorter_cols)]
            if len(line) == len(longer_cols):
                assert set(line.keys()) == longer_cols
            else:
                assert set(line.keys()) == shorter_cols
            n_lines += 1
        assert n_lines == expected_n_lines

    def test_read_variant_file_2(self, resultfile_aav2):
        
        expected_lines = {
            "reference_name":"AAV2", 
            "pos":"1485", 
            "query_name":"AAV2_N496D", 
            "var":"A1486G",
            "ref_bases":"A", 
            "query_bases":"G", 
            "aa_change":"True"
        }
        lines = [i for i in read_variant_file(resultfile_aav2)]
        assert lines == [expected_lines]

    def test_read_variant_file_3(self, resultfile_aav2_shorter):
        
        expected_lines = {
            "pos":"1485", 
            "query_name":"AAV2_N496D", 
            "ref_bases":"A", 
            "query_bases":"G", 
            "aa_change":"True"
        }
        lines = [i for i in read_variant_file(resultfile_aav2_shorter)]
        assert lines == [expected_lines]

    def test_read_variant_file_empty(self):

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            lines = [i for i in read_variant_file(temp.name)]
        assert lines == []

    def test_read_variant_file_header_only(self):

        header = "query_name\tzero_based_pos\tref_bases\tquery_bases\taa_change\n"
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            temp.write(header)
            temp.seek(0)
            lines = list(read_variant_file(temp.name))

        assert lines == []

class TestGetVariantType:

    @pytest.mark.parametrize('ref,alt,expected', 
                                [("G", "A", "sub"), 
                                 (".", "A", "ins"), 
                                 ("A", ".", "del")])
    def test_get_variant_type(self, ref, alt, expected):
        assert get_variant_type(ref, alt) == expected

class TestGetVariant:

    @pytest.mark.parametrize("write_variants", [(True, 'some_variants'), (False, 'some_variants')], indirect=['write_variants'])
    def test_get_variant(self, write_variants):
        
        _, expected_vars, temp = write_variants

        for i, line in enumerate(read_variant_file(temp.name)):
            assert get_variant(line) == expected_vars[i]

        temp.close()

class TestGetVariantsSet:

    @pytest.mark.parametrize("write_variants", [(True, 'some_variants'), (False, 'some_variants')], indirect=['write_variants'])
    def test_get_variants_set(self, write_variants):

        _, expected_vars, temp = write_variants

        expected_set = set(expected_vars)
        assert get_variants_set(temp.name) == expected_set

        temp.close()

class TestSubstituion:

    @pytest.fixture
    def substitution(self):
        rpos = 10
        rseq = 'A'
        qseq = 'T'
        changes_aa = True
        return Substitution(rpos, rseq, qseq, changes_aa)
    
    @pytest.fixture
    def substitution2(self):
        rpos = 15
        rseq = 'C'
        qseq = 'G'
        changes_aa = False
        return Substitution(rpos, rseq, qseq, changes_aa)
    
    def init_substittuion(self):
        rpos = 10
        rseq = 'A'
        qseq = 'T'
        changes_aa = True
        substitution = Substitution(rpos, rseq, qseq, changes_aa)
        return substitution

    def test_substitution(self, substitution):
        
        assert substitution.var_type == "sub"
        assert substitution.rpos == 10
        assert substitution.rseq == 'A'
        assert substitution.qseq == 'T'
        assert substitution.changes_aa == True

    def test_zero_pos(self, substitution):
        
        assert substitution.zero_pos() == 10

    def test_one_pos(self, substitution):
        
        assert substitution.one_pos() == 11

    def test_refbases(self, substitution):
        
        assert substitution.refbases() == 'A'

    def test_qbases(self, substitution):
        
        assert substitution.qbases() == 'T'

    @pytest.mark.parametrize('query_name', [None, 'query'])
    @pytest.mark.parametrize('reference_name', [None, 'reference'])
    def test_print_line(self, substitution, query_name, reference_name):
        
        if query_name is None:
            with pytest.raises(AssertionError):
                substitution.print_line(query_name, reference_name)
            return

        if reference_name is None:
            expected_line = f"{query_name}\t{substitution.zero_pos()}\t{substitution.refbases()}\t{substitution.qbases()}\t{substitution.changes_aa}\n"
            assert substitution.print_line(query_name, reference_name) == expected_line
            return
        
        expected_line = f"{reference_name}\t{substitution.zero_pos()}\t{query_name}\t{str(substitution)}\t{substitution.refbases()}\t{substitution.qbases()}\t{substitution.changes_aa}\n"
        assert substitution.print_line(query_name, reference_name) == expected_line


    def test_header(self, substitution):
        
        shorter_header = "query_name\tpos\tref_bases\tquery_bases\taa_change\n"
        longer_header = "reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n"
        
        assert substitution.header(shorter=True) == shorter_header
        assert substitution.header(shorter=False) == longer_header

    def test_var_id(self, substitution):
        
        expected_var_id = f"{substitution.zero_pos()}:{substitution.var_type}"
        
        assert substitution.var_id() == expected_var_id

    def test_str(self, substitution):
        
        expected_str = f"{substitution.rseq}{substitution.rpos+1}{substitution.qseq}"
        
        assert str(substitution) == expected_str

    def test_repr(self, substitution):
        
        expected_repr = (f"Substitution of query base '{substitution.qseq}' for reference base '{substitution.rseq}' "
                         f"at 0-based reference position {substitution.rpos}")
        
        assert repr(substitution) == expected_repr

    def test_eq(self, substitution):
        
        assert substitution == substitution

    def test_eq_not(self, substitution, substitution2):
        
        assert substitution != substitution2

    def test_hash(self, substitution):
        
        assert hash(substitution) == hash(str(substitution))

class TestInsertion:

    @pytest.fixture
    def insertion_1(self):
        last_rpos = 49
        qpos = 61
        first_query_base = 'C'
        ins = Insertion(last_rpos, qpos, first_query_base)
        return ins

    @pytest.fixture
    def insertion_3(self, insertion_1):
        last_rpos = 49
        qpos = 61
        first_query_base = 'C'
        ins = Insertion(last_rpos, qpos, first_query_base)
        ins.add_another_base('G')
        ins.add_another_base('C')
        return ins


    def test_init(self, insertion_1):
        assert insertion_1.var_type == "ins"
        assert insertion_1.last_rpos == 49
        assert insertion_1.start_qpos == 61
        assert insertion_1.end_qpos == 62
        assert insertion_1.bases == 'C'
        assert insertion_1.changes_aa == True

    def test_add_another_base(self, insertion_3):
        assert insertion_3.bases == 'CGC'
        assert insertion_3.end_qpos == 64

    @pytest.mark.parametrize('insertion', 
                             ['insertion_1', 'insertion_3'])
    def test_zero_pos(self, insertion, request):
        ins = request.getfixturevalue(insertion)
        assert ins.zero_pos() == 50

    @pytest.mark.parametrize('insertion', 
                             ['insertion_1', 'insertion_3'])
    def test_one_pos(self, insertion, request):
        ins = request.getfixturevalue(insertion)
        assert ins.one_pos() == "50_51"

    @pytest.mark.parametrize('insertion', 
                             ['insertion_1', 'insertion_3'])
    def test_refbases(self, insertion, request):
        ins = request.getfixturevalue(insertion)
        assert ins.refbases() == "."    

    @pytest.mark.parametrize('insertion,expected', 
                             [('insertion_1','C'), ('insertion_3','CGC')])
    def test_refbases(self, insertion, expected, request):
        ins = request.getfixturevalue(insertion)
        assert ins.qbases() == expected

    @pytest.mark.parametrize('insertion,expected',
                                [('insertion_1','query\t50\t.\tC\tTrue\n'), 
                                ('insertion_3','query\t50\t.\tCGC\tTrue\n')])
    def test_print_line_shorter(self, insertion, expected, request):
        ins = request.getfixturevalue(insertion)
        assert ins.print_line('query', None) == expected

    @pytest.mark.parametrize('insertion,expected',
                                [('insertion_1','reference\t50\tquery\t50_51insC\t.\tC\tTrue\n'), 
                                ('insertion_3','reference\t50\tquery\t50_51insCGC\t.\tCGC\tTrue\n')])
    def test_print_line_longer(self, insertion, expected, request):
        ins = request.getfixturevalue(insertion)
        assert ins.print_line('query', 'reference') == expected

    @pytest.mark.parametrize('shorter', [True, False])
    def test_header(self, shorter, insertion_1):
        ins = insertion_1
        shorter_header = "query_name\tpos\tref_bases\tquery_bases\taa_change\n"
        longer_header = "reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n"
        assert ins.header(shorter) == shorter_header if shorter else longer_header

    @pytest.mark.parametrize('insertion', 
                             ['insertion_1', 'insertion_3'])
    def test_var_id(self, insertion, request):
        ins = request.getfixturevalue(insertion)
        assert ins.var_id() == '50:ins'

    @pytest.mark.parametrize('insertion,expected', 
                             [('insertion_1','50_51insC'), ('insertion_3','50_51insCGC')])
    def test_str(self, insertion, expected, request):
        ins = request.getfixturevalue(insertion)
        assert ins.__str__() == expected     

    @pytest.mark.parametrize('insertion', 
                             ['insertion_1', 'insertion_3'])
    def test_repr(self, insertion, request):
        ins = request.getfixturevalue(insertion)
        expected = (f"Insertion of query bases {ins.start_qpos}:{ins.end_qpos} "
                    f" ({ins.bases}) after reference position {ins.last_rpos}")
        assert ins.__repr__() == expected      

    @pytest.mark.parametrize('insertion1,insertion2', 
                             [['insertion_1', 'insertion_1']])
    def test_eq(self, insertion1, insertion2, request):
        ins1 = request.getfixturevalue(insertion1)
        ins2 = request.getfixturevalue(insertion2)
        assert ins1 == ins2

    @pytest.mark.parametrize('insertion1,insertion2',
                                [('insertion_1', 'insertion_3'), 
                                ('insertion_1', 'insertion_3')])
    def test_eq_not(self, insertion1, insertion2, request):
        ins1 = request.getfixturevalue(insertion1)
        ins2 = request.getfixturevalue(insertion2)
        assert ins1 != ins2

    @pytest.mark.parametrize('insertion',
                                ['insertion_1', 'insertion_3'])
    def test_hash(self, insertion, request):
        ins = request.getfixturevalue(insertion)
        assert hash(ins) == hash(str(ins))

class TestDeletion:

    @pytest.fixture
    def deletion_1(self):
        rpos = 49
        last_qpos = 61
        first_ref_base = 'C'
        delete = Deletion(rpos, last_qpos, first_ref_base)
        return delete

    @pytest.fixture
    def deletion_3(self):
        rpos = 49
        last_qpos = 61
        first_ref_base = 'C'
        delete = Deletion(rpos, last_qpos, first_ref_base)
        delete.add_another_base('G')
        delete.add_another_base('C')
        return delete


    def test_init(self, deletion_1):
        assert deletion_1.var_type == "del"
        assert deletion_1.start_rpos == 49
        assert deletion_1.end_rpos == 50
        assert deletion_1.last_qpos == 61
        assert deletion_1.bases == 'C'
        assert deletion_1.changes_aa == True
        

    def test_add_another_base(self, deletion_3):
        assert deletion_3.bases == 'CGC'
        assert deletion_3.end_rpos == 52

    @pytest.mark.parametrize('deletion, expected', 
                             [('deletion_1', "49_50"), 
                              ('deletion_3', "49_52")])
    def test_zero_pos(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.zero_pos() == expected

    @pytest.mark.parametrize('deletion, expected', 
                             [('deletion_1', "50_50"), 
                              ('deletion_3', "50_52")])
    def test_one_pos(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.one_pos() == expected

    @pytest.mark.parametrize('deletion, expected', 
                             [('deletion_1', "C"), 
                              ('deletion_3', "CGC")])
    def test_refbases(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.refbases() == expected    

    @pytest.mark.parametrize('deletion', 
                             ['deletion_1', 'deletion_3'])
    def test_refbases(self, deletion, request):
        delete = request.getfixturevalue(deletion)
        assert delete.qbases() == '.'

    @pytest.mark.parametrize('deletion,expected',
                                [('deletion_1','query\t49_50\tC\t.\tTrue\n'), 
                                ('deletion_3','query\t49_52\tCGC\t.\tTrue\n')])
    def test_print_line_shorter(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.print_line('query', None) == expected

    @pytest.mark.parametrize('deletion,expected',
                                [('deletion_1','reference\t49_50\tquery\tC50del\tC\t.\tTrue\n'), 
                                 ('deletion_3','reference\t49_52\tquery\tC50_C52del\tCGC\t.\tTrue\n')])
    def test_print_line_longer(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.print_line('query', 'reference') == expected

    @pytest.mark.parametrize('shorter', [True, False])
    def test_header(self, shorter, deletion_1):
        delete = deletion_1
        shorter_header = "query_name\tpos\tref_bases\tquery_bases\taa_change\n"
        longer_header = "reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n"
        assert delete.header(shorter) == shorter_header if shorter else longer_header

    @pytest.mark.parametrize('deletion,expected',
                                [('deletion_1','49_50:del'), 
                                 ('deletion_3','49_52:del')])
    def test_var_id(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.var_id() == expected

    @pytest.mark.parametrize('deletion,expected', 
                             [('deletion_1','C50del'), 
                              ('deletion_3','C50_C52del')])
    def test_str(self, deletion, expected, request):
        delete = request.getfixturevalue(deletion)
        assert delete.__str__() == expected     

    @pytest.mark.parametrize('deletion', 
                             ['deletion_1', 'deletion_3'])
    def test_repr(self, deletion, request):
        delete = request.getfixturevalue(deletion)
        expected = (f"Deletion of reference bases {delete.start_rpos}:{delete.end_rpos} "
                    f" ({delete.bases})")    
        assert delete.__repr__() == expected      

    @pytest.mark.parametrize('deletion1,deletion2', 
                             [['deletion_1', 'deletion_1']])
    def test_eq(self, deletion1, deletion2, request):
        delete1 = request.getfixturevalue(deletion1)
        delete2 = request.getfixturevalue(deletion2)
        assert delete1 == delete2

    @pytest.mark.parametrize('deletion1,deletion2',
                                [('deletion_1', 'deletion_3'), 
                                ('deletion_1', 'deletion_3')])
    def test_eq_not(self, deletion1, deletion2, request):
        delete1 = request.getfixturevalue(deletion1)
        delete2 = request.getfixturevalue(deletion2)
        assert delete1 != delete2

    @pytest.mark.parametrize('deletion',
                                ['deletion_1', 'deletion_3'])
    def test_hash(self, deletion, request):
        delete = request.getfixturevalue(deletion)
        assert hash(delete) == hash(str(delete))