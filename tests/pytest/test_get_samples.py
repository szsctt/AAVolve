
import tempfile
import os

import pytest
import pandas as pd

from scripts.get_samples import PARENTDIR, DEFAULT_MINREPS
from scripts.get_samples import get_name, get_first_parent, get_command_options

class TestGetName:
    
    def test_get_name_simple(self):
        
        filename = "tests/data/test.fa"
        expected_result = "test"
        assert get_name(filename) == expected_result

    def test_get_name_complex(self):
        
        filename = "tests/data/test.blah.fa"
        expected_result = "test.blah"
        assert get_name(filename) == expected_result

    def test_get_name_noext(self):
        
        filename = "tests/data/test"
        expected_result = "test"
        assert get_name(filename) == expected_result

class TestFirstParent:

    def write_one_seq(self, handle, name, seq, wrapping=False):
        if not wrapping:
            handle.writelines([f">{name}\n", f"{seq}\n"])
        else:
            # split sequence into 60 character lines
            seq_split = [seq[i:i+60] for i in range(0, len(seq), 60)]
            handle.writelines([f">{name}\n"])
            for s in seq_split:
                handle.writelines([f"{s}\n"])

    @pytest.mark.parametrize("wrapping", [False, True])
    @pytest.mark.parametrize("multiple_seqs", [False, True])
    def test_first_parent(self, wrapping, multiple_seqs):

        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            
            # write sequence to temporary file
            expected_names = ['parent1', 'parent2', 'parent3']
            expected_seqs = ['ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG' * 30] * 3
            self.write_one_seq(f, expected_names[0], expected_seqs[0], wrapping=wrapping)
            
            # read written lines
            f.seek(0)
            expected_lines = f.readlines()

            # write other parents
            if multiple_seqs:
                for n, s in zip(expected_names[1:], expected_seqs[1:]):
                    self.write_one_seq(f, n, s, wrapping=wrapping)
        
            # get name and new file
            name, file = get_first_parent(f.name)

            # checks
            assert name == expected_names[0]
            assert os.path.isfile(file)
            with open(file, 'r') as p:
                lines = p.readlines()
            assert lines == expected_lines
            assert file == os.path.join(PARENTDIR, expected_names[0] + '.fa')

            # clean up
            os.remove(file)

    def test_first_parent_nofile(self):
        with pytest.raises(Exception):
            get_first_parent('nonexistentfile')

    def test_first_parent_noheader(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as f:
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
            f.seek(0)
            with pytest.raises(Exception):
                get_first_parent(f.name)

class TestGetCommandOptions:

    @pytest.mark.parametrize("sample_name", (None, 'test'))
    @pytest.mark.parametrize("parent_name", (None, 'test'))
    @pytest.mark.parametrize("reference_name", (None, 'test'))
    @pytest.mark.parametrize("seq_tech", (None, 'np', 'np-cc', 'pb'))
    @pytest.mark.parametrize("min_reps", (None, 1, 2))
    @pytest.mark.parametrize("read_file", (None, 'tests/data/reads/np-aav2.fastq'))
    @pytest.mark.parametrize("parent_file", (None, 'tests/data/references/wtAAV2.fa'))
    @pytest.mark.parametrize("reference_file", (None, 'tests/data/references/wtAAV2.fa'))
    def test_get_command_options(self, sample_name, parent_name, reference_name, seq_tech, min_reps, read_file, parent_file, reference_file):


        # create input config
        config = {}
        for k, v in zip(['sample_name', 'parent_name', 'reference_name', 'seq_tech', 'min_reps', 'read_file', 'parent_file', 'reference_file'], 
                        [sample_name, parent_name, reference_name, seq_tech, min_reps, read_file, parent_file, reference_file]):
            if v is not None:
                config[k] = v

        # we expect a fail if any of read_file, parent_file or seq_tech are missing
        if 'read_file' in config and 'parent_file' in config and 'seq_tech' in config:
            expect_fail = False
        else:
            expect_fail = True

        # inputs
        if expect_fail:
            with pytest.raises(Exception):
                samples = get_command_options(config)
            return
        else:
            samples = get_command_options(config)

        # expected defaults
        expected_config = {
            'read_file': config['read_file'], 
            'sample_name': 'np-aav2',
            'parent_file': config['parent_file'], 
            'parent_name': 'wtAAV2',
            'reference_file': os.path.join(PARENTDIR, 'AAV2.fa'), 
            'reference_name': 'AAV2',
            'seq_tech': 'np'}
        # if reference file specified breference name is derived from that by default
        if 'reference_file' in config:
            expected_config['reference_name'] = get_name(config['reference_file'])
        # update with input config
        expected_config.update(config)
        expected_samples = pd.DataFrame(expected_config, index=[0])
        
        # checks
        assert len(samples) == len(expected_samples)
        assert set(samples.columns) == set(expected_samples.columns)
        for c in samples.columns:
            assert samples[c][0] == expected_samples[c][0]

        # clean up created file
        if 'reference_file' not in config:
            os.remove(expected_config['reference_file'])


