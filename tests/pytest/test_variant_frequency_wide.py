import tempfile
import pytest

from aavolve.variant_frequency_wide import count_freqs, main

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

class TestCountFreqs:

    def test_count_freqs(self, test_wide):

        expected_freqs = {
            'var1': {'A': 0.5, 'C': 0.5},
            'var2': {'A': 0.5, 'C': 0.5},
            'var3': {'B': 0.25, 'C': 1},
            'var4': {'A': 0.5, 'C': 0.5}
        }

        freqs = count_freqs(test_wide.name, split_counts=False)

        assert freqs == expected_freqs

        test_wide.close()

    def test_count_freqs_split(self, test_wide):

        expected_freqs = {
            'var1': {'A': 0.5, 'C': 0.5},
            'var2': {'A': 0.5, 'C': 0.5},
            'var3': {'B': 0.125, 'C': 0.875},
            'var4': {'A': 0.5, 'C': 0.5}
        }

        freqs = count_freqs(test_wide.name, split_counts=True)

        assert freqs == expected_freqs

        test_wide.close()

    def test_count_freqs_header_only(self):
        
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:

            expected_freqs = {'var1' : {}, 'var2' : {}, 'var3' : {}, 'var4' : {}}

            temp.write('read_id\tvar1\tvar2\tvar3\tvar4\n')
            temp.seek(0)

            freqs = count_freqs(temp.name)

            assert freqs == expected_freqs

    def test_count_freqs_empty(self):

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:

            with pytest.raises(ValueError) as e:
                count_freqs(temp.name)
            assert str(e.value) == 'No reads found in input file'

    def test_count_freqs_no_vars(self):

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:

            temp.write('read_id\n')
            temp.seek(0)

            freqs = count_freqs(temp.name)

            assert freqs == {}

class TestMain:

    def test_main(self, monkeypatch, test_wide):

        expected_lines = [
            'variant\tparent\tfrequency\n',
            'var1\tA\t0.5\n',
            'var1\tC\t0.5\n',
            'var2\tA\t0.5\n',
            'var2\tC\t0.5\n',
            'var3\tB\t0.25\n',
            'var3\tC\t1.0\n',
            'var4\tA\t0.5\n',
            'var4\tC\t0.5\n'
        ]

        with tempfile.NamedTemporaryFile(mode='w+t') as outfile:

            monkeypatch.setattr('sys.argv', ['variant_frequency_wide.py', '-i', test_wide.name, '-o', outfile.name])

            main()

            outfile.seek(0)
            assert outfile.readlines() == expected_lines

            test_wide.close()

    def test_main_split(self, monkeypatch, test_wide):

        expected_lines = [
            'variant\tparent\tfrequency\n',
            'var1\tA\t0.5\n',
            'var1\tC\t0.5\n',
            'var2\tA\t0.5\n',
            'var2\tC\t0.5\n',
            'var3\tB\t0.125\n',
            'var3\tC\t0.875\n',
            'var4\tA\t0.5\n',
            'var4\tC\t0.5\n'
        ]

        with tempfile.NamedTemporaryFile(mode='w+t') as outfile:

            monkeypatch.setattr('sys.argv', ['variant_frequency_wide.py', '-i', test_wide.name, '-o', outfile.name, '--split-counts'])

            main()

            outfile.seek(0)
            assert outfile.readlines() == expected_lines

            test_wide.close()

    def test_main_header_only(self, monkeypatch):

        with tempfile.NamedTemporaryFile(mode='w+t') as infile, tempfile.NamedTemporaryFile(mode='w+t') as outfile:

            infile.write('read_id\tvar1\tvar2\tvar3\tvar4\n')
            infile.seek(0)

            monkeypatch.setattr('sys.argv', ['variant_frequency_wide.py', '-i', infile.name, '-o', outfile.name])

            main()

            outfile.seek(0)
            assert outfile.read() == 'variant\tparent\tfrequency\n'

    def test_main_empty(self, monkeypatch):

        with tempfile.NamedTemporaryFile(mode='w+t') as infile, tempfile.NamedTemporaryFile(mode='w+t') as outfile:

            monkeypatch.setattr('sys.argv', ['variant_frequency_wide.py', '-i', infile.name, '-o', outfile.name])

            with pytest.raises(ValueError) as e:
                main()

            assert str(e.value) == 'No reads found in input file'

    def test_main_no_vars(self, monkeypatch):

        with tempfile.NamedTemporaryFile(mode='w+t') as infile, tempfile.NamedTemporaryFile(mode='w+t') as outfile:

            infile.write('read_id\nread_1\nread_2\nread_3\n')
            infile.seek(0)

            monkeypatch.setattr('sys.argv', ['variant_frequency_wide.py', '-i', infile.name, '-o', outfile.name])

            main()

            outfile.seek(0)
            assert outfile.read() == 'variant\tparent\tfrequency\n'

