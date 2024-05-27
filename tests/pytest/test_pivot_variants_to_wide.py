import tempfile
import pytest

from aavolve.utils import Substitution, Insertion, Deletion
from aavolve.pivot_variants_to_wide import pivot_reads, get_reads, collect_read_vars, get_read_id, get_parents, closest_parent, main

@pytest.fixture
def test_group_variants():
    return [
        'tests/data/variants/test_group_variants.tsv',
        'tests/data/variants/test_group_reads.tsv',
        'tests/data/variants/test_group_parents.tsv'
    ]

@pytest.fixture
def test_group_variants_parents(test_group_variants):
    return get_parents(test_group_variants[2])

def write_header(filehandle):
    filehandle.write('reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n')
    filehandle.seek(0)

def write_read_ids(in_filename, out_handle):

    with open(in_filename, 'r') as f:
        next(f) # skip header
        read_ids = []
        for line in f:
            rid = line.split('\t')[2]
            if rid not in read_ids:
                read_ids.append(rid)

    out_handle.write('\n'.join(read_ids))
    out_handle.seek(0)

class TestPivotReads:

    def test_pivot_reads(self, resultfile_aav2389_some2_variants, resultfile_aav2389_some):
        
        parents = resultfile_aav2389_some2_variants    

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_2\tAAV3b,AAV8,AAV9\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_3\tAAV3b,AAV8,AAV9\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tT\tT\tC\tA\tC\n', # read_1 dosen't match any parents at 40:sub - could potentially be alternate allele instead of NA
            'read_2\tA\tT\tT\tC\tA\tC\n', 
            'read_3\tA\tT\tT\tC\tA\tC\n'
            ]

        # pivot the reads
        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):


            pivot_reads(resultfile_aav2389_some[0], resultfile_aav2389_some[1], outfile_parents.name, outfile_seq.name, parents, False, False, 0, 0)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_remove_na(self, resultfile_aav2389_some2_variants, resultfile_aav2389_some):

        parents = resultfile_aav2389_some2_variants
        infile, in_read_ids = resultfile_aav2389_some

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_2\tAAV3b,AAV8,AAV9\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_3\tAAV3b,AAV8,AAV9\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_2\tA\tT\tT\tC\tA\tC\n', 
            'read_3\tA\tT\tT\tC\tA\tC\n'
            ]

        # pivot the reads
        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq,
            ):

            pivot_reads(infile, in_read_ids, outfile_parents.name, outfile_seq.name, parents, True, False, 0, 0)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_no_parents(self, resultfile_aav2389_some):

        parents = {}
        infile, in_read_ids = resultfile_aav2389_some

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
        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq,
              ):

            pivot_reads(infile, in_read_ids, outfile_parents.name, outfile_seq.name, parents, False, False, 0, 0)

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
              tempfile.NamedTemporaryFile('w+t') as infile_read_ids,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            
            # write header to infile_reads
            write_header(infile_reads)

            write_read_ids(infile_reads.name, infile_read_ids)

            pivot_reads(infile_reads.name, infile_read_ids.name, outfile_parents.name, outfile_seq.name, parents, False, False, 0, 0)

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
              tempfile.NamedTemporaryFile('w+t') as infile_read_ids,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            # write header to infile_reads
            write_header(infile_reads)

            pivot_reads(infile_reads.name, infile_read_ids.name, outfile_parents.name, outfile_seq.name, parents, False, False, 0, 0)

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_group(self, test_group_variants, test_group_variants_parents):

        reads, read_ids, _ = test_group_variants

        expected_parents = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV2\n',
            'r2\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV2\n',
            'r3\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r4\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tAAV2\n',
            'r6\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r7\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r8\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\n',
        ]
        expected_seq = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tA\tA\tA\tA\tA\tA\tC\n',
            'r2\tA\tA\tA\tA\tA\tA\tC\n',
            'r3\tT\tT\tT\tT\tT\tT\tC\n',
            'r4\tT\tT\tT\tT\tT\tT\tC\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tC\n',
            'r6\tT\tT\tT\tT\tT\tT\tC\n',
            'r7\tT\tT\tT\tT\tT\tT\tC\n',
            'r8\tC\tC\tC\tC\tC\tC\tC\n',
        ]

        # pivot the reads
        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents,
                tempfile.NamedTemporaryFile('w+t') as outfile_seq
                ):
    
                pivot_reads(reads, read_ids, outfile_parents.name, outfile_seq.name, test_group_variants_parents, False, True, 1, 0.2)
    
                outfile_parents.seek(0), outfile_seq.seek(0)
                result_parents = outfile_parents.readlines()
                result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_pivot_reads_group_dist5(self, test_group_variants, test_group_variants_parents):

        reads, read_ids, _ = test_group_variants

        expected_parents = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\n',
            'r2\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\n',
            'r3\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r4\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r6\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r7\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r8\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\n',
        ]
        expected_seq = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tA\tA\tA\tA\tA\tA\tA\n',
            'r2\tA\tA\tA\tA\tA\tA\tA\n',
            'r3\tT\tT\tT\tT\tT\tT\tA\n',
            'r4\tT\tT\tT\tT\tT\tT\tA\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r6\tT\tT\tT\tT\tT\tT\tA\n',
            'r7\tT\tT\tT\tT\tT\tT\tA\n',
            'r8\tC\tC\tC\tC\tC\tC\tC\n',
        ]

        # pivot the reads
        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents,
                tempfile.NamedTemporaryFile('w+t') as outfile_seq
                ):
    
                pivot_reads(reads, read_ids, outfile_parents.name, outfile_seq.name, test_group_variants_parents, False, True, 5, 0.3)
    
                outfile_parents.seek(0), outfile_seq.seek(0)
                result_parents = outfile_parents.readlines()
                result_seq = outfile_seq.readlines()
        
        assert result_parents == expected_parents
        assert result_seq == expected_seq

class TestClosestParent:
    
    def test_closest_parent_one(self):

        group = ['40:sub', '41:sub', '42:sub', '43:sub', '44:sub', '45:sub']

        group_vars = [
            Substitution(40, 'A', 'C', True),
            Substitution(41, 'T', 'C', True),
            Substitution(42, 'G', 'C', True),
            Substitution(43, 'T', 'C', True),
            Substitution(44, 'C', 'T', False),
            Substitution(45, 'A', 'C', False),
        ]

        parents = {
            'AAV2': {},
            'AAV3b': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'C', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            },
            'AAV8': {
                '40:sub': Substitution(40, 'A', 'G', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'A', True),
                '43:sub': Substitution(43, 'T', 'A', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            }
        }

        expected_result = {            
            'AAV3b': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'C', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            },}
        
        res = closest_parent(group, group_vars, parents, 0)
        assert res == expected_result

    def test_closest_parent_split(self):

        group = ['40:sub', '41:sub', '42:sub', '43:sub', '44:sub', '45:sub']

        group_vars = [
            Substitution(40, 'A', 'C', True),
            Substitution(41, 'T', 'C', True),
            Substitution(42, 'G', 'A', True),
            Substitution(43, 'T', 'G', True),
            Substitution(44, 'C', 'T', False),
            Substitution(45, 'A', 'C', False),
        ]

        parents = {
            'AAV2': {},
            'AAV3b': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'C', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            },
            'AAV8': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'A', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            }
        }

        expected_result = {            
            'AAV3b': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'C', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            },
            'AAV8': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'A', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            }
            }

        res = closest_parent(group, group_vars, parents, 0.4)
        assert res == expected_result

    def test_closest_parent_na(self):

        group = ['40:sub', '41:sub', '42:sub', '43:sub', '44:sub', '45:sub']

        group_vars = [
            Substitution(40, 'A', 'C', True),
            Substitution(41, 'T', 'C', True),
            Substitution(42, 'G', 'A', True),
            Substitution(43, 'T', 'G', True),
            Substitution(44, 'C', 'T', False),
            Substitution(45, 'A', 'C', False),
        ]

        parents = {
            'AAV2': {},
            'AAV3b': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'C', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            },
            'AAV8': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'A', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            }
        }

        expected_result = {'NA': 'NA'}

        res = closest_parent(group, group_vars, parents, 0)
        assert res == expected_result

    def test_closest_parent_wt(self):

        group = ['40:sub', '41:sub', '42:sub', '43:sub', '44:sub', '45:sub']

        group_vars = []

        parents = {
            'AAV2': {},
            'AAV3b': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'C', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            },
            'AAV8': {
                '40:sub': Substitution(40, 'A', 'C', True),
                '41:sub': Substitution(41, 'T', 'C', True),
                '42:sub': Substitution(42, 'G', 'C', True),
                '43:sub': Substitution(43, 'T', 'A', True),
                '44:sub': Substitution(44, 'C', 'T', False),
                '45:sub': Substitution(45, 'A', 'C', False),
            }
        }

        expected_result = {'AAV2': {}}

        res = closest_parent(group, group_vars, parents, 0)
        assert res == expected_result

class TestCollectReadVars:

    def test_collect_read_vars(self):

        expected_vars = [
            ('read_1', {'40:sub': Substitution(40, 'G', 'C', True),
                        '41:sub': Substitution(41, 'T', 'C', True),
                        '44:sub': Substitution(44, 'C', 'T', False),
                        }),
            ('read_3', {'40:sub': Substitution(40, 'A', 'C', True),
                        '41:sub': Substitution(41, 'T', 'C', True),
                        }),
            ('read_4', {'40:sub': Substitution(40, 'A', 'C', True),
                        '41:sub': Substitution(41, 'T', 'C', True),
                        }),
            ('read_7', {'40:sub': Substitution(40, 'A', 'C', True),
                        '41:sub': Substitution(41, 'T', 'C', True),
                        }),
            ('read_8', {'40:sub': Substitution(40, 'G', 'C', True),
                        '41:sub': Substitution(41, 'T', 'C', True),
                        }),
            ('read_9', {'40:sub': Substitution(40, 'G', 'C', True),
                        '41:sub': Substitution(41, 'T', 'C', True),
                        }),
        ]
        with tempfile.NamedTemporaryFile('w+t') as in_vars:
            in_vars.write(expected_vars[0][1]['40:sub'].header(shorter=False))
            for read_id, rvars in expected_vars:
                for var in rvars.values():
                    in_vars.write(var.print_line(query_name=read_id, ref_name="par"))
            in_vars.seek(0)

            vars = [i for i in collect_read_vars(in_vars.name)]        

        assert vars == expected_vars

    def test_collect_read_vars_no_vars(self):

        expected_vars = []

        with tempfile.NamedTemporaryFile('w+t') as f:
            # write header to f
            write_header(f)
            f.seek(0)
            vars = [i for i in collect_read_vars(f.name)]
        
        assert vars == expected_vars
 
class TestMakeVarGroups:
    
    def test_make_var_groups(self):

        pass

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
        
        var_file = resultfile_aav2389_some[0]
        with tempfile.NamedTemporaryFile('w+t') as f:
            write_read_ids(var_file, f)

            vars = [i for i in get_reads(var_file, f.name)]

        assert len(vars) == len(expected_vars)
        assert [i[0] == j[0] and i[1]==j[1] for i, j in zip(vars, expected_vars)]

    def test_get_reads_no_reads(self):

        expected_vars = []

        with (tempfile.NamedTemporaryFile('w+t') as f,
              tempfile.NamedTemporaryFile('w+t') as f2
              ):
            # write header to f
            write_header(f)
            f.seek(0)
            vars = [i for i in get_reads(f.name, f2.name)]
        
        assert vars == expected_vars

    def test_more_readids_than_variants(self):

        vars = [
            ('read_1', {'40:sub': Substitution(40, 'G', 'C', True)}),
            ('read_3', {'40:sub': Substitution(40, 'A', 'C', True)}),
            ('read_4', {'40:sub': Substitution(40, 'A', 'C', True)}),
            ('read_7', {'40:sub': Substitution(40, 'A', 'C', True)}),
            ('read_8', {'40:sub': Substitution(40, 'G', 'C', True)}),
            ('read_9', {'40:sub': Substitution(40, 'G', 'C', True)}),
        ]
        expected_reads = [f'read_{i}' for i in range(0, 13)]
        expected_vars = vars.copy()
        expected_vars.insert(0, ('read_0', {}) )
        expected_vars.insert(2, ('read_2', {}))
        expected_vars.insert(5, ('read_5', {}))
        expected_vars.insert(6, ('read_6', {}))
        expected_vars.append(('read_10', {}))
        expected_vars.append(('read_11', {}))
        expected_vars.append(('read_12', {}))

        # write variants and read ids to file
        with (tempfile.NamedTemporaryFile('w+t') as in_vars,
              tempfile.NamedTemporaryFile('w+t') as in_read_ids,
              ):
            
            # write variants
            in_vars.write(vars[0][1]['40:sub'].header(shorter=False))
            for read_id, vars in vars:
                for var in vars.values():
                    in_vars.write(var.print_line(query_name=read_id, ref_name="par"))
            in_vars.seek(0)

            # write read ids
            in_read_ids.write('\n'.join(expected_reads))
            in_read_ids.seek(0)

            # run get_reads
            vars = [i for i in get_reads(in_vars.name, in_read_ids.name)]

            # compare
            assert vars == expected_vars
       
class TestGetReadId:

    def test_get_read_id(self):

        # write read ids
        expected_read_ids = ['read_1', 'read_2', 'read_3', 'read_4', 'read_5', 'read_6', 'read_7', 'read_8', 'read_9', 'read_10']
        with tempfile.NamedTemporaryFile('w+t') as f:
            f.write('\n'.join(expected_read_ids))
            f.seek(0)

            # get read ids
            read_ids = [i for i in get_read_id(f.name)]
        
        # check
        assert read_ids == expected_read_ids

    def test_get_read_id_no_reads(self):

        expected_read_ids = []

        # empty file
        with tempfile.NamedTemporaryFile('w+t') as f:
            read_ids = [i for i in get_read_id(f.name)]
        
        # check
        assert read_ids == expected_read_ids

class TestMain:

    def test_main(self, resultfile_aav2389_some2, resultfile_aav2389_some):
    
        infile, in_read_ids = resultfile_aav2389_some
        parents = resultfile_aav2389_some2

        expected_parents = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_2\tAAV3b,AAV8,AAV9\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n',
            'read_3\tAAV3b,AAV8,AAV9\tAAV3b,AAV9\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV8\tAAV2,AAV3b,AAV9\tAAV3b,AAV8\n'
        ]
        expected_seq = [
            'read_id\t40:sub\t44:sub\t45:sub\t46:sub\t50:sub\t53:sub\n', 
            'read_1\tNA\tT\tT\tC\tA\tC\n', # read_1 dosen't match any parents at 40:sub - could potentially be alternat allele instead of NA
            'read_2\tA\tT\tT\tC\tA\tC\n', 
            'read_3\tA\tT\tT\tC\tA\tC\n'
            ]

        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq
              ):

            main(['-i', infile, '-r', in_read_ids, '-p', parents,'-o', outfile_parents.name, '-O', outfile_seq.name])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_no_parents(self, resultfile_aav2389_some):
    
        infile, in_read_ids = resultfile_aav2389_some

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

            main(['-i', infile,  '-r', in_read_ids, '-p', infile_parents.name,
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
              tempfile.NamedTemporaryFile('w+t') as infile_read_ids,
              tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq):
            # write header to infile_reads
            write_header(infile_reads)

            write_read_ids(infile_reads.name, infile_read_ids)

            main(['-i', infile_reads.name, '-r', infile_read_ids.name,
                  '-p', parents, '-o', outfile_parents.name,
                  '-O', outfile_seq.name])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_grouping(self, test_group_variants):

        reads, read_ids, parents = test_group_variants

        expected_parents = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV2\n',
            'r2\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV2\n',
            'r3\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r4\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tAAV2\n',
            'r6\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r7\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV2\n',
            'r8\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\n',
        ]
        expected_seq = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tA\tA\tA\tA\tA\tA\tC\n',
            'r2\tA\tA\tA\tA\tA\tA\tC\n',
            'r3\tT\tT\tT\tT\tT\tT\tC\n',
            'r4\tT\tT\tT\tT\tT\tT\tC\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tC\n',
            'r6\tT\tT\tT\tT\tT\tT\tC\n',
            'r7\tT\tT\tT\tT\tT\tT\tC\n',
            'r8\tC\tC\tC\tC\tC\tC\tC\n',
        ]

        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
              tempfile.NamedTemporaryFile('w+t') as outfile_seq
              ):

            main(['-i', reads, '-r', read_ids, '-p', parents, '-o', outfile_parents.name, '-O', outfile_seq.name, '-g'])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_grouping_combine(self, test_group_variants):

        reads, read_ids, parents = test_group_variants

        expected_parents = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\n',
            'r2\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r3\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r4\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r6\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r7\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r8\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\n',
        ]
        expected_seq = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tA\tA\tA\tA\tA\tA\tA\n',
            'r2\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r3\tT\tT\tT\tT\tT\tT\tA\n',
            'r4\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r6\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r7\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r8\tC\tC\tC\tC\tC\tC\tC\n',
        ]

        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
            tempfile.NamedTemporaryFile('w+t') as outfile_seq
            ):

            main(['-i', reads, '-r', read_ids, '-p', parents, '-o', outfile_parents.name, '-O', outfile_seq.name, '-g', '--group-dist', '5'])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_grouping_combine_na(self, test_group_variants):

        reads, read_ids, parents = test_group_variants

        expected_parents = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\n',
            'r3\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r8\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\n',
        ]
        expected_seq = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tA\tA\tA\tA\tA\tA\tA\n',
            'r3\tT\tT\tT\tT\tT\tT\tA\n',
            'r8\tC\tC\tC\tC\tC\tC\tC\n',
        ]

        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
            tempfile.NamedTemporaryFile('w+t') as outfile_seq
            ):

            main(['-i', reads, '-r', read_ids, '-p', parents, '-o', outfile_parents.name, '-O', outfile_seq.name, '-g', '--group-dist', '5', '--remove-na'])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq

    def test_main_grouping_combine_freq(self, test_group_variants):

        reads, read_ids, parents = test_group_variants

        expected_parents = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\n',
            'r2\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\tAAV3b\n',
            'r3\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r4\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r6\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r7\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\tAAV8\n',
            'r8\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\tAAV2\n',
        ]
        expected_seq = [
            'read_id\t40:sub\t41:sub\t42:sub\t43:sub\t44:sub\t45:sub\t50:sub\n',
            'r1\tA\tA\tA\tA\tA\tA\tA\n',
            'r2\tA\tA\tA\tA\tA\tA\tA\n',
            'r3\tT\tT\tT\tT\tT\tT\tA\n',
            'r4\tT\tT\tT\tT\tT\tT\tA\n',
            'r5\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n',
            'r6\tT\tT\tT\tT\tT\tT\tA\n',
            'r7\tT\tT\tT\tT\tT\tT\tA\n',
            'r8\tC\tC\tC\tC\tC\tC\tC\n',
        ]

        with (tempfile.NamedTemporaryFile('w+t') as outfile_parents, 
            tempfile.NamedTemporaryFile('w+t') as outfile_seq
            ):

            main(['-i', reads, '-r', read_ids, '-p', parents, '-o', outfile_parents.name, '-O', outfile_seq.name, '-g', '--group-dist', '5', '--max-distance-frac', '0.5'])

            outfile_parents.seek(0), outfile_seq.seek(0)
            result_parents = outfile_parents.readlines()
            result_seq = outfile_seq.readlines()

        assert result_parents == expected_parents
        assert result_seq == expected_seq