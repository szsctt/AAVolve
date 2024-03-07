import tempfile
import pytest
from Bio import SeqIO
from Bio import Seq

from scripts.apply_variants import read_reference, apply_variants, create_temp_seq_file, main
from scripts.utils import use_open

TOY_EXCLUDED = {'read1_matches_ref', 
                'read3_sub_A6G_insbeforeStart', 
                'read13_sub_C1G_A98_A100del',
                'read14_C1_T5del', 
                'read15_C1_T5del_A98_A100del', 
                'read16_1_insaaa_C1_T5del_A98_A100del'}


def read_fa(filename):
    d = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    d = {k: str(v.seq).upper() for k, v in d.items()}
    return d

def read_variants(filename):
     
    with use_open(filename,  'rt') as f:
        header = f.readline().strip().split('\t')

        for row in f:
            yield dict(zip(header, row.strip().split('\t')))


class TestImportReadReference:

    def test_read_reference(self, aav2_ref_file):

        expected_ref = 'ATGGCTGCCGATGGTTATCTTCCAGATTGGCTCGAGGACACTCTCTCTGAAGGAATAAGACAGTGGTGGAAGCTCAAACCTGGCCCACCACCACCAAAGCCCGCAGAGCGGCATAAGGACGACAGCAGGGGTCTTGTGCTTCCTGGGTACAAGTACCTCGGACCCTTCAACGGACTCGACAAGGGAGAGCCGGTCAACGAGGCAGACGCCGCGGCCCTCGAGCACGACAAAGCCTACGACCGGCAGCTCGACAGCGGAGACAACCCGTACCTCAAGTACAACCACGCCGACGCGGAGTTTCAGGAGCGCCTTAAAGAAGATACGTCTTTTGGGGGCAACCTCGGACGAGCAGTCTTCCAGGCGAAAAAGAGGGTTCTTGAACCTCTGGGCCTGGTTGAGGAACCTGTTAAGACGGCTCCGGGAAAAAAGAGGCCGGTAGAGCACTCTCCTGTGGAGCCAGACTCCTCCTCGGGAACCGGAAAGGCGGGCCAGCAGCCTGCAAGAAAAAGATTGAATTTTGGTCAGACTGGAGACGCAGACTCAGTACCTGACCCCCAGCCTCTCGGACAGCCACCAGCAGCCCCCTCTGGTCTGGGAACTAATACGATGGCTACAGGCAGTGGCGCACCAATGGCAGACAATAACGAGGGCGCCGACGGAGTGGGTAATTCCTCGGGAAATTGGCATTGCGATTCCACATGGATGGGCGACAGAGTCATCACCACCAGCACCCGAACCTGGGCCCTGCCCACCTACAACAACCACCTCTACAAACAAATTTCCAGCCAATCAGGAGCCTCGAACGACAATCACTACTTTGGCTACAGCACCCCTTGGGGGTATTTTGACTTCAACAGATTCCACTGCCACTTTTCACCACGTGACTGGCAAAGACTCATCAACAACAACTGGGGATTCCGACCCAAGAGACTCAACTTCAAGCTCTTTAACATTCAAGTCAAAGAGGTCACGCAGAATGACGGTACGACGACGATTGCCAATAACCTTACCAGCACGGTTCAGGTGTTTACTGACTCGGAGTACCAGCTCCCGTACGTCCTCGGCTCGGCGCATCAAGGATGCCTCCCGCCGTTCCCAGCAGACGTCTTCATGGTGCCACAGTATGGATACCTCACCCTGAACAACGGGAGTCAGGCAGTAGGACGCTCTTCATTTTACTGCCTGGAGTACTTTCCTTCTCAGATGCTGCGTACCGGAAACAACTTTACCTTCAGCTACACTTTTGAGGACGTTCCTTTCCACAGCAGCTACGCTCACAGCCAGAGTCTGGACCGTCTCATGAATCCTCTCATCGACCAGTACCTGTATTACTTGAGCAGAACAAACACTCCAAGTGGAACCACCACGCAGTCAAGGCTTCAGTTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAAAGACATCTGCGGATAACAACAACAGTGAATACTCGTGGACTGGAGCTACCAAGTACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGCCCGGCCATGGCAAGCCACAAGGACGATGAAGAAAAGTTTTTTCCTCAGAGCGGGGTTCTCATCTTTGGGAAGCAAGGCTCAGAGAAAACAAATGTGGACATTGAAAAGGTCATGATTACAGACGAAGAGGAAATCAGGACAACCAATCCCGTGGCTACGGAGCAGTATGGTTCTGTATCTACCAACCTCCAGAGAGGCAACAGACAAGCAGCTACCGCAGATGTCAACACACAAGGCGTTCTTCCAGGCATGGTCTGGCAGGACAGAGATGTGTACCTTCAGGGGCCCATCTGGGCAAAGATTCCACACACGGACGGACATTTTCACCCCTCTCCCCTCATGGGTGGATTCGGACTTAAACACCCTCCTCCACAGATTCTCATCAAGAACACCCCGGTACCTGCGAATCCTTCGACCACCTTCAGTGCGGCAAAGTTTGCTTCCTTCATCACACAGTACTCCACGGGACAGGTCAGCGTGGAGATCGAGTGGGAGCTGCAGAAGGAAAACAGCAAACGCTGGAATCCCGAAATTCAGTACACTTCCAACTACAACAAGTCTGTTAATGTGGACTTTACTGTGGACACTAATGGCGTGTATTCAGAGCCTCGCCCCATTGGCACCAGATACCTGACTCGTAATCTGTAA'
        ref = read_reference(aav2_ref_file)

        assert ref == expected_ref

    def test_read_reference_multiple(self, aav2389_ref_file): 
        """
        Should error for more than one reference
        """

        with pytest.raises(AssertionError):
            read_reference(aav2389_ref_file)

class TestApplyVariants:

    def test_apply_variants_1(self):
        """
        test with made-up example

        num: 0123456789.012345678901234567890123456789
        ref: ATGCATGCAT.GCATGCATGCATGCATGCATGCATGCATGC
        mut: ATGCACGCATTGCATG.ATGCATGCATGCATGCATGCATGC
                  ^    ^     ^
        """

        row = {'5:sub': 'C', '10:ins': 'T', '15_16:del': '.'}
  
        ref =          'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC'
        expected_seq = 'ATGCACGCATTGCATGATGCATGCATGCATGCATGCATGC'

        count, seq = apply_variants(ref, row)

        assert count == 1
        assert seq == expected_seq

    def test_apply_variants_2(self):
        """
        test with made-up example where substitution overlaps
        with deletion, which can happen when subuttion is
        reference allele and deletion is not or vice versa 
                       1         2         3         
        num: 0123456789012345678901234567890123456789
        ref: ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
        mut: ATGCACGCATGCATG.....ATGCATGCATGCATGCATGC
                  ^         ^^^^^
        """

        row = {'5:sub': 'C', '15_20:del': '.', '17:sub': 'T'}
  
        ref =          'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC'
        expected_seq = 'ATGCACGCATGCATGATGCATGCATGCATGCATGC'

        count, seq = apply_variants(ref, row)

        assert count == 1
        assert seq == expected_seq

    def test_apply_variants_round_trip_1(self, aav3_pivoted_seq, aav23_ref_file):
        """
        Get variants of AAV3 relative to AAV2, then apply them to AAV2 to get back to AAV3
        """

        refs = read_fa(aav23_ref_file)

        with open(aav3_pivoted_seq,  'rt') as f:
            header = f.readline().strip().split('\t')
            row = f.readline().strip().split('\t')
            variants = dict(zip(header, row))

        count, seq = apply_variants(refs['AAV2'], variants)

        assert seq == refs['AAV3']
        assert count == 1
        
    def test_apply_variants_round_trip_2(self, toy_ref_file, toy_reads_file, toy_pivoted_seq):
        """
        Get variants of toy reads relative to reference, then apply to get back to reads
        """

        ref = read_fa(toy_ref_file)
        expected = read_fa(toy_reads_file)

        for row in read_variants(toy_pivoted_seq):

            read_id = row['read_id']
            expected_seq = expected[read_id]

            # insertions before start of reference are not reproduced
            if read_id in TOY_EXCLUDED:
                continue
            
            count, seq = apply_variants(ref['ref'], row)
                
            assert seq == expected_seq
            assert count == 1


class TestCreateTempSeqFile:
    
    def test_create_temp_seq_file_aav3(self, aav3_pivoted_seq, aav3_pivoted_parents, resultfile_aav23):
        """
        Test that create_temp_seq_file creates a file with the expected sequence
        """

        # get expected file
        with use_open(aav3_pivoted_seq, 'rt') as f:
            expected = f.readlines()

        # run fucntion and read result
        res = create_temp_seq_file(aav3_pivoted_parents, resultfile_aav23)
        result = res.readlines()

        # check
        assert result == expected

        # clean up
        res.close()

    def test_create_seq_file_toy(self, toy_pivoted_seq, toy_pivoted_parents, resultfile_toy):
        """
        Test that create_temp_seq_file creates a file with the expected sequence
        """

        # get expected file
        with use_open(toy_pivoted_seq, 'rt') as f:
            expected = f.readlines()

        # run function and read result
        res = create_temp_seq_file(toy_pivoted_parents, resultfile_toy)
        result = res.readlines()

        # check
        assert result == expected

        # clean up
        res.close()

        
class TestMain:

    @pytest.mark.parametrize('translate', (True, False))
    @pytest.mark.parametrize('from_parent_names', (True, False))
    @pytest.mark.parametrize('fasta_output', (True, False))
    def test_main_toy_seq(self, toy_ref_file,toy_pivoted_seq, toy_pivoted_parents, resultfile_toy, toy_reads_file, 
                          monkeypatch, translate, from_parent_names, fasta_output):
        """
        Test main function
        """
        
        # not all reads are included in the results 
        with open(toy_pivoted_seq, 'rt') as f:
            expected_reads = f.readlines()
            expected_reads = [i.split('\t')[0] for i in expected_reads[1:] if i.split('\t')[0]]

        # get sequences of reads
        reads = read_fa(toy_reads_file)
        expected_seqs = [i for k, i in reads.items() if k in expected_reads]
        # translate if expected
        if translate:
            # pad to length that is multiple of three
            to_add = [3 - len(i) % 3 for i in expected_seqs]
            to_add = [i if i < 3 else 0 for i in to_add]
            expected_seqs = [i + 'N' * j for i, j in zip(expected_seqs, to_add)]
            expected_seqs = [Seq.translate(i) for i in expected_seqs]
        # insertions before the start of the read aren't reproduced, so read 3 is the same as read 2
        expected_seqs[1] = expected_seqs[0]
        
        # reproduce expected output
        if fasta_output:
            expected = []
            for i, seq in enumerate(expected_seqs):
                expected.append(f'>seq{i}_1\n')
                expected.append(f'{seq}\n')
        else:
            expected = ['count\tsequence\n'] + [f'1\t{i}\n' for i in expected_seqs]
        
        # run main
        with tempfile.NamedTemporaryFile(mode='w+t') as f:

            # set sys.argv with monkeypatch
            if from_parent_names:
                args = ['main', '-v', toy_pivoted_parents, '-p', resultfile_toy, '-r', toy_ref_file, '-o', f.name]
            else:
                args = ['main', '-v', toy_pivoted_seq, '-r', toy_ref_file, '-o', f.name]
            if translate:
                args.append('--translate')
            if fasta_output:
                args.append('--fasta')
            monkeypatch.setattr('sys.argv', args)

            main()
        
            # read results
            f.seek(0)
            result = f.readlines()

        # check
        assert result == expected
  

class TestMain:

    @pytest.mark.parametrize('translate', (True, False))
    @pytest.mark.parametrize('from_parent_names', (True, False))
    @pytest.mark.parametrize('fasta_output', (True, False))
    def test_main_aav3(self,aav3_pivoted_seq, aav3_pivoted_parents, resultfile_aav23, aav23_ref_file, aav2_ref_file,
                          monkeypatch, translate, from_parent_names, fasta_output):
        """
        Test main function
        """
        
        # get sequence of AAV2 and AAV3
        refs = read_fa(aav23_ref_file)

        expected_seq = refs['AAV3']
        # translate if expected
        if translate:
            expected_seq = Seq.translate(expected_seq)
        
        # reproduce expected output
        if fasta_output:
            expected = [f'>seq0_1\n', f'{expected_seq}\n']
        else:
            expected = ['count\tsequence\n', f'1\t{expected_seq}\n'] 
        
        # run main
        with tempfile.NamedTemporaryFile(mode='w+t') as f:

            # set sys.argv with monkeypatch
            if from_parent_names:
                args = ['main', '-v', aav3_pivoted_parents, '-p', resultfile_aav23, '-r', aav2_ref_file, '-o', f.name]
            else:
                args = ['main', '-v', aav3_pivoted_seq, '-r', aav2_ref_file, '-o', f.name]
            if translate:
                args.append('--translate')
            if fasta_output:
                args.append('--fasta')
            monkeypatch.setattr('sys.argv', args)

            main()
        
            # read results
            f.seek(0)
            result = f.readlines()

        # check
        assert result == expected
  