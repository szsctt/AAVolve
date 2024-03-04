import pytest
from Bio import SeqIO

from scripts.apply_variants import read_reference, apply_variants, create_temp_seq_file, main
from scripts.utils import use_open


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
            if read_id in {'read3_sub_A6G_insbeforeStart'}:
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
        res = create_temp_seq_file(resultfile_aav23, aav3_pivoted_parents)
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
        res = create_temp_seq_file(resultfile_toy, toy_pivoted_parents)
        result = res.readlines()

        # check
        assert result == expected

        # clean up
        res.close()

        


        