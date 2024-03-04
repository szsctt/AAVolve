import tempfile
import pytest

from scripts.utils import Substitution, Insertion, Deletion
from scripts.pivot_variants_to_wide import pivot_reads, get_reads, main

def write_header(filehandle):
    filehandle.write('reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n')
    filehandle.seek(0)


class TestPivotReads:

    def test_pivot_reads(self, resultfile_aav2389_some2_variants, resultfile_aav2389_some):
        
        parents = resultfile_aav2389_some2_variants
        infile = resultfile_aav2389_some

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_2\tAAV3b,AAV9,AAV8\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_3\tAAV3b,AAV9,AAV8\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tT\tT\tC\tA\tC\n', # read_1 dosen't match any parents at 40:sub - could potentially be alternat allele instead of NA
            'read_2\tA\tT\tT\tC\tA\tC\n', 
            'read_3\tA\tT\tT\tC\tA\tC\n'
            ]

        # pivot the reads
        with tempfile.NamedTemporaryFile('w+t') as outfile_parents, tempfile.NamedTemporaryFile('w+t') as outfile_seq:
            pivot_reads(infile, outfile_parents.name, outfile_seq.name, parents, False)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_remove_na(self, resultfile_aav2389_some2_variants, resultfile_aav2389_some):

        parents = resultfile_aav2389_some2_variants
        infile = resultfile_aav2389_some

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_2\tAAV3b,AAV9,AAV8\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_3\tAAV3b,AAV9,AAV8\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_2\tA\tT\tT\tC\tA\tC\n', 
            'read_3\tA\tT\tT\tC\tA\tC\n'
            ]

        # pivot the reads
        with tempfile.NamedTemporaryFile('w+t') as outfile_parents, tempfile.NamedTemporaryFile('w+t') as outfile_seq:
            pivot_reads(infile, outfile_parents.name, outfile_seq.name, parents, True)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_no_parents(self, resultfile_aav2389_some):

        parents = {}
        infile = resultfile_aav2389_some

        expected_parents = [
            'read_id\n', 
            'read_1\n',
            'read_2\n',
            'read_3\n'
        ]
        expected_seq = [
            'read_id\n', 
            'read_1\n', 
            'read_2\n', 
            'read_3\n'
            ]

        # pivot the reads
        with tempfile.NamedTemporaryFile('w+t') as outfile_parents, tempfile.NamedTemporaryFile('w+t') as outfile_seq:
            pivot_reads(infile, outfile_parents.name, outfile_seq.name, parents, False)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_no_reads(self, resultfile_aav2389_some2_variants):

        parents = resultfile_aav2389_some2_variants

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
        ]

        # pivot the reads
        with (tempfile.NamedTemporaryFile('w+t') as infile_reads,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            
            # write header to infile_reads
            write_header(infile_reads)

            pivot_reads(infile_reads.name, outfile_parents.name, outfile_seq.name, parents, False)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq


    def test_pivot_reads_no_parents_or_reads(self):

        parents = {}

        expected_parents = [
            'read_id\n'
        ]
        expected_seq = [
            'read_id\n', 
        ]

        # pivot the reads
        with (tempfile.NamedTemporaryFile('w+t') as infile_reads,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            # write header to infile_reads
            write_header(infile_reads)

            pivot_reads(infile_reads.name, outfile_parents.name, outfile_seq.name, parents, False)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

class TestGetReads:

    def test_get_reads(self, resultfile_aav2389_some):

        expected_vars = [
                ('read_1', {'40:sub': Substitution(40, 'G', 'C', True),
                            '41:sub': Substitution(41, 'T', 'C', True),
                            '44:sub': Substitution(44, 'C', 'T', False),
                            '53:sub': Substitution(53, 'A', 'C', False),
                            '55:sub': Substitution(55, 'A', 'C', False),
                            }),
                ('read_2', {'40:sub': Substitution(40, 'C', 'A', True),
                            '41:sub': Substitution(41, 'T', 'C', True),
                            '44:sub': Substitution(44, 'C', 'T', False),
                            '53:sub': Substitution(53, 'A', 'C', False),
                            '55:sub': Substitution(55, 'A', 'C', False),
                            }),
                ('read_3', {'40:sub': Substitution(40, 'C', 'A', True),
                            '41:sub': Substitution(41, 'T', 'C', True),
                            '44:sub': Substitution(44, 'C', 'T', False),
                            '53:sub': Substitution(53, 'A', 'C', False),
                }),                 
        ]
        
        var_file = resultfile_aav2389_some

        vars = [i for i in get_reads(var_file)]

        assert len(vars) == len(expected_vars)
        assert [i[0] == j[0] and i[1]==j[1] for i, j in zip(vars, expected_vars)]

    def test_get_reads_no_reads(self):

        expected_vars = []

        with tempfile.NamedTemporaryFile('w+t') as f:
            # write header to f
            write_header(f)
            f.seek(0)
            vars = [i for i in get_reads(f.name)]
        
        assert vars == expected_vars


class TestMain:

    def test_main(self, resultfile_aav2389_some2, resultfile_aav2389_some):
    
        infile = resultfile_aav2389_some
        parents = resultfile_aav2389_some2

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_2\tAAV3b,AAV9,AAV8\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_3\tAAV3b,AAV9,AAV8\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tT\tT\tC\tA\tC\n', # read_1 dosen't match any parents at 40:sub - could potentially be alternat allele instead of NA
            'read_2\tA\tT\tT\tC\tA\tC\n', 
            'read_3\tA\tT\tT\tC\tA\tC\n'
            ]

        with tempfile.NamedTemporaryFile('w+t') as outfile_parents, tempfile.NamedTemporaryFile('w+t') as outfile_seq:

            main(['-i', infile, '-p', parents,'-o', outfile_parents.name, '-O', outfile_seq.name])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_no_parents(self, resultfile_aav2389_some):
    
        infile = resultfile_aav2389_some

        expected_parents = [
            'read_id\n', 
            'read_1\n',
            'read_2\n',
            'read_3\n'
        ]
        expected_seq = [
            'read_id\n', 
            'read_1\n', 
            'read_2\n', 
            'read_3\n'
            ]

        with (tempfile.NamedTemporaryFile('w+t') as infile_parents,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            
            # write header for parents file
            write_header(infile_parents)


            main(['-i', infile, '-p', infile_parents.name,
                   '-o', outfile_parents.name, '-O', outfile_seq.name])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_no_reads(self, resultfile_aav2389_some2):
    
        parents = resultfile_aav2389_some2

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
        ]

        with (tempfile.NamedTemporaryFile('w+t') as infile_reads,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            # write header to infile_reads
            write_header(infile_reads)

            main(['-i', infile_reads.name, '-p', parents, '-o', outfile_parents.name,
                  '-O', outfile_seq.name])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq