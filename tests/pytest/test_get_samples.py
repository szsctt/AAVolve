
import tempfile
import os

import pytest
import pandas as pd

from scripts.get_samples import PARENTDIR, DEFAULT_MINREPS, REQUIRED_COLUMNS, SEQ_TECHS
from scripts.get_samples import get_name, get_first_parent, get_command_options, check_data, get_samples

class Fixtures:

    @pytest.fixture
    def config(self):
        return {
            'sample_name': 'sample1',
            'parent_name': 'aavs',
            'reference_name': 'AAV2',
            'seq_tech': 'np',
            'read_file': 'tests/data/reads/np-aav2.fastq',
            'parent_file': 'tests/data/references/AAV2_AAV3.fa',
            'reference_file': 'tests/data/references/wtAAV2.fa'
        }
    
    @pytest.fixture 
    def sample_df(self, config):

        return pd.DataFrame([config])

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
        # if reference file specified reference name is derived from that by default
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

class TestCheckData(Fixtures):
    
    def test_check_data(self, sample_df):
        """
        Minimal test to check that no exception is raised
        """

        # should not raise
        check_data(sample_df)

    @pytest.mark.parametrize("col", REQUIRED_COLUMNS)
    def test_check_data_missing_column(self, sample_df, col):
        """
        Check that an exception is raised if a required column is missing
        """

        # remove required column
        sample_df = sample_df.drop(columns=[col])

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == f'Missing required column: {col}'

    @pytest.mark.parametrize("col", REQUIRED_COLUMNS)
    def test_check_data_missing_value(self, sample_df, col):
        """
        Check that an exception is raised if a required value is missing
        """

        # remove required value
        sample_df[col] = None

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == f'Missing value in column: {col}'

    @pytest.mark.parametrize("col", ("read_file", "parent_file", "reference_file"))
    def test_check_data_missing_file(self, sample_df, col):
        """
        Check that an exception is raised if a required file does not exist
        """

        # remove required file
        sample_df[col] = 'nonexistentfile'

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)
        
        # check error message
        types = {'read_file': 'Sample', 'parent_file': 'Parent', 'reference_file': 'Reference'}
        assert error.value.args[0] == f'{types[col]} file does not exist: {sample_df[col][0]}'

    def test_check_data_duplicate_samples(self, sample_df):
        """
        Check that an exception is raised if sample names are not unique
        """

        # duplicate sample
        sample_df = pd.concat([sample_df, sample_df])

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == "Sample names (column 'sample_name') must be unique"

    def test_check_data_ref_files(self, sample_df):
        """
        Check that an exception is raised if different reference_file used for the same reference_name
        """

        # same reference name but different reference files
        sample_df = pd.concat([sample_df, sample_df])
        sample_df['sample_name'] = ['sample1', 'sample2']
        sample_df['reference_file'] = ['tests/data/references/wtAAV2.fa', 'tests/data/references/AAV2_AAV3.fa']

        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == "Each reference name (column 'reference_name') must always correspond to the same reference file (column 'reference_file')"

    def test_check_data_ref_names(self, sample_df):
        """
        Check that an exception is raised if different reference_name used for the same reference_file
        """

        # same reference file but different reference names
        sample_df = pd.concat([sample_df, sample_df])
        sample_df['sample_name'] = ['sample1', 'sample2']
        sample_df['reference_name'] = ['AAV2', 'AAV3']

        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == "Each reference name (column 'reference_name') must always correspond to the same reference file (column 'reference_file')"

    def test_check_data_parent_files(self, sample_df):
        """
        Check that an exception is raised if different parent_file used for the same parent_name
        """
    
        # same parent name but different parent files
        sample_df = pd.concat([sample_df, sample_df])
        sample_df['sample_name'] = ['sample1', 'sample2']
        sample_df['parent_file'] = ['tests/data/references/wtAAV2.fa', 'tests/data/references/AAV2_AAV3.fa']

        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == "Each parent name (column 'parent_name') must always correspond to the same parent file (column 'parent_file')"

    def test_check_data_parent_names(self, sample_df):
        """
        Check that an exception is raised if different parent_name used for the same parent_file
        """

        # same parent file but different parent names
        sample_df = pd.concat([sample_df, sample_df])
        sample_df['sample_name'] = ['sample1', 'sample2']
        sample_df['parent_name'] = ['aavs', 'aav2']

        with pytest.raises(Exception) as error:
            check_data(sample_df)

        assert error.value.args[0] =="Each parent name (column 'parent_name') must always correspond to the same parent file (column 'parent_file')"

    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_check_data_seq_tech_exists(self, sample_df, seq_tech):
        """
        Check that no exception is raised if seq_tech is one of the allowed values
        """

        # change seq tech
        sample_df['seq_tech'] = seq_tech

        # should not raise
        check_data(sample_df)

    def test_check_data_seq_tech_not_exists(self, sample_df):
        """
        Check that an exception is raised if seq_tech is not one of the allowed values
        """
            
        # change seq tech
        sample_df['seq_tech'] = 'nonexistent'

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)

        assert error.value.args[0] == "Sequencing technology (column 'seq_tech') must be one of 'np', 'np-cc', 'pb', 'pb-hifi'"

    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_check_data_min_reps_added(self, sample_df, seq_tech):
        """
        Check min_reps column is added with correct default if not present
        """

        # add min reps
        sample_df['seq_tech'] = seq_tech

        # should not raise
        check_data(sample_df)

        assert 'min_reps' in sample_df.columns
        if seq_tech == 'np-cc':
            assert sample_df['min_reps'][0] == DEFAULT_MINREPS
        else:
            assert sample_df['min_reps'][0] is None

    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_check_data_min_reps_negative(self, sample_df, seq_tech):
        """
        Check that an exception is raised if min_reps is negative
        """

        # change min reps
        sample_df['seq_tech'] = seq_tech
        sample_df['min_reps'] = -1

        # if seq_tech is np-cc, should raise
        if seq_tech == 'np-cc':
            with pytest.raises(Exception) as error:
                check_data(sample_df)
            assert error.value.args[0] == "Minimum reps (column 'min_reps') must be at least 0"
        # otherwise set to null
        else:
            check_data(sample_df)
            assert all(sample_df['min_reps'].isnull())

    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_check_data_min_reps_null(self, sample_df, seq_tech):
        """
        Check that min_reps is set to default if null and seq_tech is np-cc, otherwise set to null
        """

        # change min reps
        sample_df['seq_tech'] = seq_tech
        sample_df['min_reps'] = None

        # should not raise
        check_data(sample_df)

        if seq_tech == 'np-cc':
            assert sample_df['min_reps'][0] == DEFAULT_MINREPS
        else:
            assert sample_df['min_reps'][0] is None

class TestGetSamples(Fixtures):

    def test_get_samples_from_command_line(self, config):
        """
        Check that get_samples returns the correct dataframe when called from the command line through --config
        argument
        """

        samples = get_samples(config)

        expected_samples = pd.DataFrame([config])
        expected_samples['min_reps'] = None 

        # check data frames are equivalent - columns might be in different order
        assert set(samples.columns) == set(expected_samples.columns)
        for col in samples.columns:
            assert [i == j for i, j in zip(samples[col], expected_samples[col])]


    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_get_samples_from_file(self, config, seq_tech):
        """
        Check that get_samples returns the correct dataframe when called from a file
        """

        # write config to file
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.csv') as f:
        
            config['seq_tech'] = seq_tech
            expected_samples = pd.DataFrame([config])
            expected_samples.to_csv(f.name, index=False)

            expected_samples['min_reps'] = DEFAULT_MINREPS if seq_tech == 'np-cc' else None

            # pass in both file and config
            config['samples'] = f.name
            config['sample_name'] = "sample2" # so we can tell if the df is created from config or file

            # read config from file
            samples = get_samples(config)

            # check data frames are equivalent - columns might be in different order
            assert set(samples.columns) == set(expected_samples.columns)
            for col in samples.columns:
                assert [i == j for i, j in zip(samples[col], expected_samples[col])]

    def test_get_samples_from_file_missing_file(self, config):

        # pass in both file and config
        config['samples'] = 'nonexistentfile'
        config['sample_name'] = "sample2"

        # should raise
        with pytest.raises(Exception) as error:
            get_samples(config)
