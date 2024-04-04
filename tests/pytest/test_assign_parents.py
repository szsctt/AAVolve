import tempfile
import pytest

from aavolve.assign_parents import get_parents, main

@pytest.fixture
def test_wide():
    temp = tempfile.NamedTemporaryFile(mode='w+t')
    temp.write('read_id\tvar1\tvar2\tvar3\tvar4\n')
    temp.write('read1\tA\tA\tB,C\tA\n')
    temp.write('read2\tA\tA\tC\tA\n')
    temp.write('read3\tC\tC\tC\tC\n')
    temp.write('read4\tC\tC\tC\tC\n')
    temp.seek(0)
    return temp

class TestAssignParents:
    
    def test_get_parents_1(self):
        
        row = ["read_1", "A,B,C", "A,C", "B,C", "A"]
        expected = ["read_1", "C", "C", "C", "A"]
        assert get_parents(row) == expected

    def test_get_parents_no_variants(self):

        row = ["read_1"]
        expected = ["read_1"]
        assert get_parents(row) == expected

    def test_get_parents_with_na(self):
        
        # NA should be ignored
        row = ["read_1", "A,B,C", "NA", "B,C", "A"]
        expected = ["read_1", "B,C", "B,C", "B,C", "A"]
        assert get_parents(row) == expected

    def test_get_parents_all_na(self):
            
            # NA should be ignored
            row = ["read_1", "NA", "NA", "NA", "NA"]
            expected = ["read_1", "NA", "NA", "NA", "NA"]
            assert get_parents(row) == expected

class TestMain:

    def test_main_1(self, monkeypatch):

        with tempfile.NamedTemporaryFile(mode="w+") as infile, tempfile.NamedTemporaryFile(mode="w+") as outfile:
            
            infile.write("read\tvariant_1\tvariant_2\tvariant_3\n")
            infile.write("read_1\tA,B,C\tA,C\tB,C\tA\n")
            infile.write("read_2\tA,B\tA,C\tA,B,C\tA\n")
            infile.write("read_3\tA,B,C\tA,C\tB,C\tA,C\n")
            infile.seek(0)

            # set sys.argv
            monkeypatch.setattr("sys.argv", ["assign_parents.py", '-i', infile.name, '-o', outfile.name])

            main()

            outfile.seek(0)
            assert outfile.read() == "read\tvariant_1\tvariant_2\tvariant_3\nread_1\tC\tC\tC\tA\nread_2\tA\tA\tA\tA\nread_3\tC\tC\tC\tC\n"

    def test_main_empty(self, monkeypatch):

        with tempfile.NamedTemporaryFile(mode="w+") as infile, tempfile.NamedTemporaryFile(mode="w+") as outfile:
            
            infile.write("read\tvariant_1\tvariant_2\tvariant_3\n")
            infile.seek(0)

            # set sys.argv
            monkeypatch.setattr("sys.argv", ["assign_parents.py", '-i', infile.name, '-o', outfile.name])

            main()

            outfile.seek(0)
            assert outfile.read() == "read\tvariant_1\tvariant_2\tvariant_3\n"