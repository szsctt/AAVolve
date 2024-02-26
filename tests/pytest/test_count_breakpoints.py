import tempfile
import pytest

from scripts.count_breakpoints import count_breakpoints, find_breakpoints, get_var_set, main

@pytest.fixture
def test_wide():
    temp = tempfile.NamedTemporaryFile(mode='w+t')
    temp.write('read_id\t1:sub\t2:sub\t5:sub\t60_71:del\t81:sub\t100_101:ins\n')
    temp.write('read1\tA\tA\tB,C\tA\tB\tB\n')
    temp.write('read2\tA\tA\tC\tA\tB,C\tC\n')
    temp.write('read3\tC\tC\tC\tC\tC\tC\n')
    temp.write('read4\tC\tC\tC\tC\tB,C\tB,A\n')
    temp.write('read5\tA,B\tB,C\tC,D\tA,C,D\tC\tC\n')
    temp.seek(0)
    return temp



class TestFindBreakpoints:

    def test_find_breakpoints(self):

        # test data
        line = {'read_id': 'read1', 'var1': 'A', 'var2': 'A', 'var3': 'B', 'var4': 'B', 'var5': 'A'}
        sorted_cols = ['var1', 'var2', 'var3', 'var4', 'var5']

        expected_breaks = ['var3', 'var5']

        # run function
        breaks = find_breakpoints(line, sorted_cols)

        # assert
        assert breaks == expected_breaks

    def test_find_breakpoints_with_ambig(self):


        # test data
        line = {'read_id': 'read1', 'var1': 'A', 'var2': 'A', 'var3': 'B', 'var4': 'A,B', 'var5': 'A'}
        sorted_cols = ['var1', 'var2', 'var3', 'var4', 'var5']

        expected_breaks = ['var3', 'var5']

        # run function
        breaks = find_breakpoints(line, sorted_cols)

        # assert
        assert breaks == expected_breaks

class TestCountBreakpoints:

    def test_count_breakpoints(self, test_wide):

        # test data
        expected_lines = [
            'read_id\t1:sub\t2:sub\t5:sub\t60_71:del\t81:sub\t100_101:ins\n',
            'read1\t0\t0\t1\t1\t1\t0\n',
            'read2\t0\t0\t1\t1\t1\t0\n',
            'read3\t0\t0\t0\t0\t0\t0\n',
            'read4\t0\t0\t0\t0\t0\t1\n',
            'read5\t0\t0\t1\t0\t0\t0\n'
        ]
        expected_n_break_reads = {'read1': 3, 'read2': 3, 'read3': 0, 'read4': 1, 'read5': 1}
        expected_n_break_loc = {'1:sub': 0, '2:sub': 0, '5:sub': 3, '60_71:del': 2, '81:sub': 2, '100_101:ins': 1}

        # run function
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            n_break_reads, n_break_loc = count_breakpoints(test_wide.name, temp.name)

            temp.seek(0)
            out_lines = temp.readlines()

        # assert
        assert out_lines == expected_lines
        assert n_break_reads == expected_n_break_reads
        assert n_break_loc == expected_n_break_loc

        test_wide.close()


class TestGetVarSet:

    def test_get_var_set(self):
            
        # test data
        vars = 'A,B,C'

        # run function
        var_set = get_var_set(vars)

        # assert
        assert var_set == {'A', 'B', 'C'}

    def test_get_var_set_2(self):

        # test data
        vars = 'A'

        # run function
        var_set = get_var_set(vars)

        # assert
        assert var_set == {'A'}

    def test_get_var_set_3(self):

        vars = 'A,B,C,NA'

        with pytest.raises(AssertionError):
            get_var_set(vars)

class TestMain:

    def test_main(self, test_wide, monkeypatch):

        # test data
        expected_lines = [
            'read_id\t1:sub\t2:sub\t5:sub\t60_71:del\t81:sub\t100_101:ins\n',
            'read1\t0\t0\t1\t1\t1\t0\n',
            'read2\t0\t0\t1\t1\t1\t0\n',
            'read3\t0\t0\t0\t0\t0\t0\n',
            'read4\t0\t0\t0\t0\t0\t1\n',
            'read5\t0\t0\t1\t0\t0\t0\n'
        ]
        expected_summary1 = [
            'read_id\tbreakpoints\n',
            'read1\t3\n',
            'read2\t3\n',
            'read3\t0\n',
            'read4\t1\n',
            'read5\t1\n'
        ]
        expected_summary2 = [
            'location\tbreakpoints\n',
            '1:sub\t0\n',
            '2:sub\t0\n',
            '5:sub\t3\n',
            '60_71:del\t2\n',
            '81:sub\t2\n',
            '100_101:ins\t1\n'
       ]

        # run function
        with (tempfile.NamedTemporaryFile(mode='w+t') as out_all,
              tempfile.NamedTemporaryFile(mode='w+t') as out_summary1,
              tempfile.NamedTemporaryFile(mode='w+t') as out_summary2):
            
            # set sys.argv
            monkeypatch.setattr("sys.argv", ["count_breakpoints.py", '-i', test_wide.name, '-o', out_all.name, '-s1', out_summary1.name, '-s2', out_summary2.name])
            
            main()

            out_all.seek(0), out_summary1.seek(0), out_summary2.seek(0)
            out_lines = out_all.readlines()
            out_summary1_lines = out_summary1.readlines()
            out_summary2_lines = out_summary2.readlines()

        # assert
        assert out_lines == expected_lines
        assert out_summary1_lines == expected_summary1
        assert out_summary2_lines == expected_summary2

        test_wide.close()
