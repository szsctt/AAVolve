import tempfile
import random
import subprocess

import pytest
import numpy as np

from scripts.distance_matrix import alignment, dmat, read_input, main

class TestReadInput:

    def test_read_input(self):

        expected_seqs = ['ACDEFGHIKLMNPQRSTVWW', 'ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVYY']

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # write header
            f.write('seq_id\tseq\n')
            
            # write sequences
            for i, seq in enumerate(expected_seqs):
                f.write(f'seq_{i}\t{seq}\n')
            f.flush()

            # read sequences
            seqs = read_input(f.name, 3, 'random')
        
        # check
        assert len(seqs) == 3
        assert all([i in seqs for i in expected_seqs])

    def test_read_input_count_one_column(self):

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # write header
            f.write('seq\n')
            
            # write sequences
            f.write('AGTC\n')
            f.flush()

            # read sequences
            seqs = read_input(f.name, 1, 'random')
        
        # check
        assert len(seqs) == 1
        assert seqs[0] == 'AGTC'

    @pytest.mark.parametrize('selection', ['random', 'first', 'last'])
    @pytest.mark.parametrize('max_seqs', [5, 10])
    def test_read_input_sample(self, selection, max_seqs):
            
        # generate 100 random sequences of characters of length 10
        sequences = [''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=10)) for _ in range(100)]

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # write header
            f.write('seq_id\tseq\n')
            
            # write sequences
            for i, seq in enumerate(sequences):
                f.write(f'seq_{i}\t{seq}\n')
            f.flush()

            # read sequences
            seqs = read_input(f.name, max_seqs, selection)

        # check
        assert len(seqs) == max_seqs
        assert all([i in sequences for i in seqs])

        if selection == 'first':
            assert seqs == sequences[:max_seqs]
        elif selection == 'last':
            assert seqs == sequences[-max_seqs:]
        
            
    def test_read_input_invalid_selection(self):

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # read sequences
            with pytest.raises(ValueError):
                read_input(f.name, 1, 'invalid')

    def test_read_input_invalid_max_seqs(self):

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # read sequences
            with pytest.raises(ValueError):
                read_input(f.name, -1, 'random')

    def test_read_input_header_only(self):

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # write header only
            f.write('seq_id\tseq\n')
            f.flush()

            # read sequences
            seqs = read_input(f.name, 3, 'random')
        
        # check
        assert len(seqs) == 0

    def test_read_input_empty(self):

        with tempfile.NamedTemporaryFile(mode='w') as f:

            # read sequences
            with pytest.raises(ValueError):
                read_input(f.name, 3, 'random')

class TestAlignment:

    def seq_comparison(self, seq1, seq2):

        return all([i==j for i,j in zip(seq1, seq2)])

    def test_alignment(self):

        # create sequences
        seqs = ['ACDEFGHIILMNPQRSTVYY', 'ACDEFGHIKLMNNQRSTVWY', 'ACDEEFGHIKLMNPQRSTVWY']
        expected_aln = ['ACD-EFGHIILMNPQRSTVYY', 'ACD-EFGHIKLMNNQRSTVWY', 'ACDEEFGHIKLMNPQRSTVWY']

        # do alignment
        aln = alignment(seqs)

        # check
        assert len(aln) == 3
        assert all([self.seq_comparison(i, j) for i, j  in zip(aln, expected_aln)])

    def test_alignment_wrong_chars(self):

        seqs = ['AGTAATGAGA', 'AGTAACJIQAGA', 'AGTAATGAGA']

        with pytest.raises(subprocess.CalledProcessError):
            alignment(seqs)

    def test_alignment_empty(self):

        # create sequences
        seqs = []

        # do alignment
        aln = alignment(seqs)

        # check
        assert len(aln) == 0

class TestDmat:

    def test_dmat_nt(self):

        seqs = ['AAAA', 'AAAC', 'AAAT']

        result = dmat(seqs, 'identity')

        assert result.shape == (3, 3)
        assert np.isclose(
            result, np.array(
            [[0.  , 0.25, 0.25],
             [0.25, 0.  , 0.25],
             [0.25, 0.25, 0.  ]]
             )
        ).all()

    def test_dmat_aa_identity(self):

        seqs = ['ACDEFGHIKLMNPQRSTVYY', 'ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVYY']

        result = dmat(seqs, 'identity')

        assert result.shape == (3, 3)
        assert np.isclose(
            result, np.array(
            [[0.  , 0.05, 0.  ],
             [0.05, 0.  , 0.05],
             [0.  , 0.05, 0.  ]]
             )
        ).all()

    def test_dmat_aa_blosum62(self):

        seqs = ['ACDEFGHIKLMNPQRSTVYY', 'ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVYY']

        result = dmat(seqs, 'blosum62')

        assert result.shape == (3, 3)
        assert np.isclose(
            result, np.array(
            [[0.  ,       0.07758621, 0.  ],
             [0.07758621, 0.  ,       0.07758621],
             [0.  ,       0.07758621, 0.  ]]
             )
        ).all()

class TestMain:

    @pytest.mark.parametrize('do_plot', (True, False))
    @pytest.mark.parametrize('metric,seqs,expected', (
            ('identity', ['AAAA', 'AAAC', 'AAAT'], np.array([[0., 0.25, 0.25], [0.25, 0.  , 0.25], [0.25, 0.25, 0.  ]])),
            ('identity', ['ACDEFGHIKLMNPQRSTVYY', 'ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVYY'], np.array([[0., 0.05, 0.  ], [0.05, 0. , 0.05], [0., 0.05, 0. ]])),
            ('blosum62', ['ACDEFGHIKLMNPQRSTVYY', 'ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVYY'], np.array([[0., 0.07758621, 0.], [0.07758621, 0., 0.07758621], [0., 0.07758621, 0.]])),
    ))
    def test_main(self, metric, seqs, expected, do_plot, monkeypatch):

        with tempfile.NamedTemporaryFile(mode='w') as infile, tempfile.NamedTemporaryFile(mode='w') as outfile, tempfile.NamedTemporaryFile(mode='w', suffix='.png') as plotfile:

            # write header
            infile.write('seq_id\tseq\n')
            
            # write sequences
            for i, seq in enumerate(seqs):
                infile.write(f'seq_{i}\t{seq}\n')
            infile.flush()

            # set sys.argv with monkeypatch
            args = ['distance_matrix.py', '-i', infile.name, '-o', outfile.name, '-d', metric, '-x', 'first']
            if do_plot:
                args.append('-p')
                args.append(plotfile.name)
            monkeypatch.setattr('sys.argv', args)

            # run main
            main()

            # check matrix result
            result = np.loadtxt(outfile.name, delimiter='\t')
            assert np.isclose(result, expected).all()

            # check plot
            with open(plotfile.name, 'rb') as f:
                contents = f.read()
                if do_plot:
                    # check first bytes are 137 80 78 71 13 10 26 10
                    # http://www.libpng.org/pub/png/spec/1.2/PNG-Structure.html
                    assert contents[:8] == b'\x89PNG\r\n\x1a\n'
                else:
                    assert contents == b''

