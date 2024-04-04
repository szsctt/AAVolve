from io import StringIO
import pytest

from aavolve.sum_counts import main

class TestMain:

    def test_main(self, monkeypatch, capsys):
    
        # set sys.stdin
        monkeypatch.setattr('sys.stdin', StringIO("1\tA\n1\tA\n1\tA\n1\tB\n1\tB\n1\tC\n"))

        # run main
        main()

        # check output
        captured = capsys.readouterr()
        assert captured.out == "3\tA\n2\tB\n1\tC\n"

    def test_main_empty(self, monkeypatch, capsys):

        # set sys.stdin
        monkeypatch.setattr('sys.stdin', StringIO(""))

        # run main
        main()

        # check output
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_main_single(self, monkeypatch, capsys):

        # set sys.stdin
        monkeypatch.setattr('sys.stdin', StringIO("1\tA\n"))

        # run main
        main()

        # check output
        captured = capsys.readouterr()
        assert captured.out == "1\tA\n"

    def test_main_with_header(self, monkeypatch):

        # set sys.stdin
        monkeypatch.setattr('sys.stdin', StringIO("count\tsequence\n1\tA\n1\tA\n1\tA\n1\tB\n1\tB\n1\tC\n"))

        # run main
        with pytest.raises(ValueError):
            main()

    def test_main_unsorted(self, monkeypatch):

        # set sys.stdin
        monkeypatch.setattr('sys.stdin', StringIO("1\tA\n1\tB\n2\tA\n1\tB\n1\tB\n1\tC\n"))

        # run main
        with pytest.raises(ValueError):
            main()