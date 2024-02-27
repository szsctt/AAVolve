
import tempfile
import os

import pytest
import pandas as pd

from scripts.get_samples import (PARENTDIR, DEFAULT_MINREPS, DEFAULT_FREQ, 
                                 DEFAULT_INCLUDE_NON, REQUIRED_COLUMNS, SEQ_TECHS)
from scripts.get_samples import get_name, get_first_parent, get_command_options, check_data, get_samples


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
    @pytest.mark.parametrize("splint_file", (None, 'tests/data/references/splint.fa'))
    def test_get_command_options(self, sample_name, parent_name, reference_name, seq_tech, min_reps, read_file, parent_file, reference_file, splint_file):


        # create input config
        config = {}
        for k, v in zip(['sample_name', 'parent_name', 'reference_name', 'seq_tech', 'min_reps', 'read_file', 'parent_file', 'reference_file', 'splint_file'], 
                        [sample_name, parent_name, reference_name, seq_tech, min_reps, read_file, parent_file, reference_file, splint_file]):
            if v is not None:
                config[k] = v

        # we expect a fail if any of read_file, parent_file or seq_tech are missing
        expect_fail = False
        if 'read_file' not in config or 'parent_file' not in config or 'seq_tech' not in config:
            expect_fail = True
        if seq_tech == 'np-cc' and 'splint_file' not in config:
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

class TestCheckData:
    
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

    def test_check_data_sample_parent_names(self, sample_df):
        """
        Check that an exception is raised if sample names are the same as parent names
        """

        # same parent name as sample name
        sample_df['parent_name'] = sample_df['sample_name']

        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == "Sample names (column 'sample_name') must be different from parent names (column 'parent_name')"

    def test_check_data_sample_reference_names(self, sample_df):
        """
        Check that an exception is raised if sample names are the same as reference names
        """

        # same reference name as sample name
        sample_df['reference_name'] = sample_df['sample_name']

        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == "Sample names (column 'sample_name') must be different from reference names (column 'reference_name')"

    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_check_data_seq_tech_exists(self, sample_df, seq_tech):
        """
        Check that no exception is raised if seq_tech is one of the allowed values
        """

        # change seq tech
        sample_df['seq_tech'] = seq_tech

        if seq_tech == 'np-cc':
            sample_df['splint_file'] = 'tests/data/references/splint.fa'

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
        
        if seq_tech == 'np-cc':
            sample_df['splint_file'] = ['tests/data/references/splint.fa']

        # should not raise
        check_data(sample_df)

        assert 'min_reps' in sample_df.columns
        if seq_tech == 'np-cc':
            assert sample_df['min_reps'][0] == DEFAULT_MINREPS
        else:
            assert sample_df['min_reps'][0] is None

    @pytest.mark.parametrize("min_reps", (-1, 'foo'))
    @pytest.mark.parametrize("seq_tech", SEQ_TECHS)
    def test_check_data_min_reps_invalid(self, sample_df, seq_tech, min_reps):
        """
        Check that an exception is raised if min_reps is negative or not a number
        """
        expected_error = f"Minimum reps (column 'min_reps') must be an integer and at least 0: found value {min_reps} in row 0"

        # change min reps
        sample_df['seq_tech'] = seq_tech
        sample_df['min_reps'] = min_reps
        if seq_tech == 'np-cc':
            sample_df['splint_file'] = ['tests/data/references/splint.fa']

        # if seq_tech is np-cc, should raise
        if seq_tech == 'np-cc':
            with pytest.raises(Exception) as error:
                check_data(sample_df)
            assert error.value.args[0] == expected_error
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
        if seq_tech == 'np-cc':
            sample_df['splint_file'] = ['tests/data/references/splint.fa']

        # should not raise
        check_data(sample_df)

        if seq_tech == 'np-cc':
            assert sample_df['min_reps'][0] == DEFAULT_MINREPS
        else:
            assert sample_df['min_reps'][0] is None

    def test_check_data_np_cc_splint(self, config):
        """
        Check that an exception is raised if splint file is not specified for np-cc
        """

        # change seq tech
        config['seq_tech'] = 'np-cc'

        # should raise
        with pytest.raises(Exception) as error:
            check_data(pd.DataFrame([config]))
        assert error.value.args[0] == "For sequencing technology 'np-cc', must specify splint file (column 'splint_file')"

    def test_check_data_np_cc_splint_exists(self, config):

        # change seq tech
        config['seq_tech'] = 'np-cc'
        config['splint_file'] = 'tests/data/references/splint.fa'

        # should not raise
        check_data(pd.DataFrame([config]))

    def test_check_data_np_cc_splint_missing_file(self, config):

        # change seq tech
        config['seq_tech'] = 'np-cc'
        config['splint_file'] = 'nonexistentfile'

        # should raise
        with pytest.raises(Exception) as error:
            check_data(pd.DataFrame([config]))
        assert error.value.args[0] == "Splint file does not exist: nonexistentfile"

    def test_check_data_non_parental_freq_added(self, sample_df):
        """
        Check that non_parental_freq is set to default if not present
        """

        check_data(sample_df)

        assert 'non_parental_freq' in sample_df.columns
        assert all(sample_df['non_parental_freq'] == DEFAULT_FREQ)

           
    @pytest.mark.parametrize("non_parental_freq", [0, 0.2, 0.5, 1])
    def test_check_data_non_parental_freq_correct(self, sample_df, non_parental_freq):
        """
        Check that an exception is raised if non_parental_freq is not between 0 and 1
        """

        # change non_parental_freq
        sample_df['non_parental_freq'] = non_parental_freq

        check_data(sample_df)

        assert 'non_parental_freq' in sample_df.columns
        assert all([i == non_parental_freq for i in sample_df['non_parental_freq']])


    @pytest.mark.parametrize("non_parental_freq", [-0.1, 1.1])
    def test_check_data_non_parental_freq_incorrect(self, sample_df, non_parental_freq):
        """
        Check that an exception is raised if non_parental_freq is not between 0 and 1
        """
        
        expected_error = f"Column 'non_parental_freq' must be numeric and between 0 and 1: found value {non_parental_freq} in row 0"

        # change non_parental_freq
        sample_df['non_parental_freq'] = non_parental_freq

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == expected_error

    def test_check_data_include_non_parental_added(self, sample_df):
        """
        Check that include_non_parental is set to default if not present
        """

        check_data(sample_df)

        assert 'include_non_parental' in sample_df.columns
        assert all(sample_df['include_non_parental'] == DEFAULT_INCLUDE_NON)
    
    @pytest.mark.parametrize("include_non_parental", [True, False])
    def test_check_data_include_non_parental_correct(self, sample_df, include_non_parental):
        """
        Check that an exception is not raised if include_non_parental not True or False
        """

        # change include_non_parental
        sample_df['include_non_parental'] = include_non_parental

        check_data(sample_df)

        assert 'include_non_parental' in sample_df.columns
        assert all([i == include_non_parental for i in sample_df['include_non_parental']])
    
    @pytest.mark.parametrize("include_non_parental", ["blah", 0.2, 1.1])
    def test_check_data_include_non_parental_incorrect(self, sample_df, include_non_parental):
        """
        Check that an exception is raised if include_non_parental not True or False
        """

        expected_error = f"Column 'include_non_parental' must be True, False or omitted: found value {include_non_parental} in row 0"

        # change include_non_parental
        sample_df['include_non_parental'] = include_non_parental

        # should raise
        with pytest.raises(Exception) as error:
            check_data(sample_df)
        assert error.value.args[0] == expected_error

class TestGetSamples:

    def test_get_samples_from_command_line(self, config):
        """
        Check that get_samples returns the correct dataframe when called from the command line through --config
        argument
        """

        samples = get_samples(config)

        expected_samples = pd.DataFrame([config])
        expected_samples['min_reps'] = None 
        expected_samples['non_parental_freq'] = DEFAULT_FREQ
        expected_samples['include_non_parental'] = DEFAULT_INCLUDE_NON

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
            if seq_tech == 'np-cc':
                config['splint_file'] = 'tests/data/references/splint.fa'
            expected_samples = pd.DataFrame([config])
            expected_samples.to_csv(f.name, index=False)

            expected_samples['min_reps'] = DEFAULT_MINREPS if seq_tech == 'np-cc' else None
            expected_samples['non_parental_freq'] = DEFAULT_FREQ
            expected_samples['include_non_parental'] = DEFAULT_INCLUDE_NON


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
