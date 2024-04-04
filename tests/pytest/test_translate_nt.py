import tempfile
import pytest

from aavolve.translate_nt import main


@pytest.fixture
def seq_counts(tmpdir):
    with tempfile.NamedTemporaryFile(mode='w+t', dir=tmpdir) as f:
        f.write('count\tsequence\n1\tATGCGT\n1\tATGCGT\n1\tATGCGT\n')
        f.seek(0)
        yield f.name, 'count\tsequence\n1\tMR\n1\tMR\n1\tMR\n'

@pytest.fixture
def seq_counts_2(tmpdir):
    with tempfile.NamedTemporaryFile(mode='w+t', dir=tmpdir) as f:
        f.write('count\tsequence\n1\tATGCGTA\n1\tATGCGT\n1\tATGCGT\n')
        f.seek(0)
        yield f.name, 'count\tsequence\n1\tMR\n1\tMR\n1\tMR\n'       

@pytest.fixture
def seq_counts_empty(tmpdir):
    with tempfile.NamedTemporaryFile(mode='w+t', dir=tmpdir) as f:
        f.write('count\tsequence\n')
        f.seek(0)
        yield f.name, 'count\tsequence\n'  

@pytest.fixture
def seq_counts_invalid(tmpdir):
    with tempfile.NamedTemporaryFile(mode='w+t', dir=tmpdir) as f:
        f.write('count\n1\t')
        f.seek(0)
        yield f.name, None


@pytest.fixture
def input_file(request):
    return request.getfixturevalue(request.param)

class TestMain:

    @pytest.mark.parametrize('input_file', ['seq_counts', 'seq_counts_2', 'seq_counts_empty', 'seq_counts_invalid'], indirect=True)
    def test_main(self, input_file, monkeypatch):

        filename, expected = input_file
        
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            monkeypatch.setattr('sys.argv', ['translate_nt.py', '-i', filename, '-o', f.name])

            if expected is None:
                with pytest.raises(AssertionError):
                    main()
                return

            main()
            f.seek(0)
            assert f.read() == expected