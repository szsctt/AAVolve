import tempfile

from scripts.remove_first_column import main

class TestMain:

    def test_main(self, monkeypatch, capsys):

        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write("a\tb\tc\n1\t2\t3\n4\t5\t6\n")
            f.seek(0)
            monkeypatch.setattr('sys.argv', ['remove_first_column.py', '-i', f.name])
            main()
            captured = capsys.readouterr()
        
        assert captured.out == "2\t3\n5\t6\n"

    def test_main_empty(self, monkeypatch, capsys):

        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write("a\n")
            f.seek(0)
            monkeypatch.setattr('sys.argv', ['remove_first_column.py', '-i', f.name])
            main()
            captured = capsys.readouterr()
        
        assert captured.out == ''

    def test_main_first_col_only(self, monkeypatch, capsys):

        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write("a\n1\n4\n")
            f.seek(0)
            monkeypatch.setattr('sys.argv', ['remove_first_column.py', '-i', f.name])
            main()
            captured = capsys.readouterr()
        
        assert captured.out == ''