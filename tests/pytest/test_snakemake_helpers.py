from types import SimpleNamespace
import pytest
import numpy as np
import pandas as pd
from aavolve.snakemake_helpers import (
    fill_parents,
    format_input_reads,
    get_column_by_parent,
    get_column_by_sample,
    get_dmat_input,
    get_linked_adapters_for_sample,
    get_parents,
    get_reads_for_align,
    get_reads_for_counting,
    get_reference,
    is_fastq,
    minimap2_params_with_default,
    trim_is_enabled,
)

class TestGetColumnBySample:
    
    def test_get_column_by_sample(self):
        df = pd.DataFrame({'sample_name': [1, 2, 3], 'sample2': [4, 5, 6]})
        wildcards = SimpleNamespace(sample=2)
        result = get_column_by_sample(wildcards, df, 'sample2')
        assert result == 5

    def test_get_column_by_sample_invalid(self):
        df = pd.DataFrame({'sample_name': [1, 2, 3], 'sample2': [4, 5, 6]})
        wildcards = SimpleNamespace(sample=4)
        with pytest.raises(KeyError):
            get_column_by_sample(wildcards, df, 'sample2')

    def test_get_column_by_sample_repeated(self):
        df = pd.DataFrame({'sample_name': [1, 2, 2], 'sample2': [4, 5, 6]})
        wildcards = SimpleNamespace(sample=2)
        with pytest.raises(AssertionError):
            get_column_by_sample(wildcards, df, 'sample2')

class TestGetColumnByParent:

    def test_get_column_by_parent_valid(self):
        wildcards = SimpleNamespace(sample='parent1')
        samples = pd.DataFrame({
            'parent_name': ['parent1', 'parent2'],
            'some_col': [42, 99]
        })
        assert get_column_by_parent(wildcards, samples, 'some_col') == 42

    def test_get_column_by_parent_not_found(self):
        wildcards = SimpleNamespace(sample='parent3')
        samples = pd.DataFrame({
            'parent_name': ['parent1', 'parent2'],
            'some_col': [42, 99]
        })
        with pytest.raises(KeyError):
            get_column_by_parent(wildcards, samples, 'some_col')

    def test_get_column_by_parent_nonunique(self):
        wildcards = SimpleNamespace(sample='parent1')
        samples = pd.DataFrame({
            'parent_name': ['parent1', 'parent1'],
            'some_col': [42, 43]
        })
        with pytest.raises(AssertionError):
            get_column_by_parent(wildcards, samples, 'some_col')

    def test_get_column_by_parent_nonunique_2(self):
        # Test with a valid parent but non-unique values in the column - this is ok
        wildcards = SimpleNamespace(sample='parent1')
        samples = pd.DataFrame({
            'parent_name': ['parent1', 'parent2'],
            'some_col': [42, 42]
        })
        assert get_column_by_parent(wildcards, samples, 'some_col') == 42

    def test_get_column_by_parent_column_missing(self):
        wildcards = SimpleNamespace(sample='parent1')
        samples = pd.DataFrame({
            'parent_name': ['parent1', 'parent2'],
            # 'some_col' missing
        })
        with pytest.raises(KeyError):
            get_column_by_parent(wildcards, samples, 'some_col')

class TestIsFastq:

    @pytest.mark.parametrize('name, res', [
        ('file.fastq', True),
        ('file.fastq.gz', True),
        ('file.fq', True),
        ('file.fq.gz', True),
        ('file.txt', False),
        ('file.txt.gz', False),
        ('file', False),
    ])
    def test_is_fastq(self, name, res):
        assert is_fastq(name) == res


class TestTrimmingHelpers:

    def _make_samples(self, **overrides):
        base = {
            'sample_name': ['sample1'],
            'parent_name': ['parent1'],
            'parent_file': ['parent1.fa'],
            'reference_file': ['ref.fa'],
            'read_file': ['reads.fastq.gz'],
            'seq_tech': ['np'],
            'min_reps': [np.nan],
            'trim': [False],
            'adapter_5': ['AAA'],
            'adapter_3': ['TTT'],
        }
        for key, value in overrides.items():
            base[key] = value
        return pd.DataFrame(base)

    def test_trim_is_enabled_false_by_default(self):
        samples = self._make_samples()
        wildcards = SimpleNamespace(sample='sample1')
        assert trim_is_enabled(wildcards, samples) is False

    def test_trim_is_enabled_false_for_parent(self):
        samples = self._make_samples(trim=[True])
        wildcards = SimpleNamespace(sample='parent1')
        assert trim_is_enabled(wildcards, samples) is False

    def test_trim_is_enabled_raises_on_non_boolean(self):
        samples = self._make_samples(trim=['YES'])
        wildcards = SimpleNamespace(sample='sample1')
        with pytest.raises(ValueError):
            trim_is_enabled(wildcards, samples)

    def test_get_reads_for_align_returns_original_when_trim_disabled(self):
        samples = self._make_samples()
        wildcards = SimpleNamespace(sample='sample1')
        reads = get_reads_for_align(wildcards, samples)
        assert reads == 'reads.fastq.gz'

    def test_get_reads_for_align_returns_trimmed_when_enabled(self):
        samples = self._make_samples(trim=[True])
        wildcards = SimpleNamespace(sample='sample1')
        reads = get_reads_for_align(wildcards, samples)
        assert reads == 'out/trimmed/sample1.trimmed.gz'

    def test_get_linked_adapters_for_sample_disabled(self):
        samples = self._make_samples()
        wildcards = SimpleNamespace(sample='sample1')
        assert get_linked_adapters_for_sample(wildcards, samples) is None

    def test_get_linked_adapters_for_sample_enabled(self):
        samples = self._make_samples(trim=[True], adapter_5=[' AAA '], adapter_3=['TTT  '])
        wildcards = SimpleNamespace(sample='sample1')
        assert get_linked_adapters_for_sample(wildcards, samples) == 'AAA...TTT'

    def test_get_linked_adapters_for_sample_missing_sequence(self):
        samples = self._make_samples(trim=[True], adapter_5=[''], adapter_3=['TTT'])
        wildcards = SimpleNamespace(sample='sample1')
        with pytest.raises(ValueError):
            get_linked_adapters_for_sample(wildcards, samples)


class TestGetReference:

    @pytest.mark.parametrize("sample, exp", [
        ('parent1', 'ref1.fa'),
        ('parent2', 'ref2.fa'),
        ('sample1', 'ref1.fa'),
        ('sample2', 'ref2.fa'),
        ('sample3', 'ref1.fa'),
        ('not a sample', 'error'),
    ])
    def test_get_reference(self, sample, exp):

        wildcards = SimpleNamespace(sample=sample)
        samples = pd.DataFrame({'parent_name': ['parent1', 'parent2', 'parent1'],
                                'sample_name': ['sample1', 'sample2', 'sample3'],
                                'reference_file': ['ref1.fa', 'ref2.fa', 'ref1.fa']})
        if exp == 'error':
            with pytest.raises(ValueError):
                get_reference(wildcards, samples)
        else:
            result = get_reference(wildcards, samples)
            assert result == exp

class TestMinimap2ParamsWithDefault:

    @pytest.mark.parametrize("params, exp", [
        ('', '-x map-hifi -B 1.5 --end-bonus 5'),
        ('-k 15', '-k 15 -x map-hifi -B 1.5 --end-bonus 5'),
        ('-w 10', '-w 10 -x map-hifi -B 1.5 --end-bonus 5'),
        ('-B 1.5', '-B 1.5 -x map-hifi --end-bonus 5'),
        ('-k 15 -w 10 -B 1.5', '-k 15 -w 10 -B 1.5 -x map-hifi --end-bonus 5'),
        ('--end-bonus 10', '--end-bonus 10 -x map-hifi -B 1.5'),
        ('-x map-ont', '-x map-ont -B 1.5 --end-bonus 5'),
        ('--preset map-ont', '--preset map-ont -B 1.5 --end-bonus 5'),
        ('-x map-ont -B 2', '-x map-ont -B 2 --end-bonus 5'),
        ('-x map-ont --end-bonus 7', '-x map-ont --end-bonus 7 -B 1.5'),
        ('-x map-ont -B 2 --end-bonus 7', '-x map-ont -B 2 --end-bonus 7'),
    ])
    def test_minimap2_params_with_default(self, params, exp):
        wildcards = SimpleNamespace(sample='sample1')
        samples = pd.DataFrame({'parent_name': ['parent1', 'parent2', 'parent1'],
                                'sample_name': ['sample1', 'sample2', 'sample3'],
                                'reference_file': ['ref1.fa', 'ref2.fa', 'ref1.fa'],
                                'minimap2_params': [params, params, params]})
        assert minimap2_params_with_default(wildcards, samples) == exp

    def test_missing_minimap2_params_column(self):
        wildcards = SimpleNamespace(sample='sample1')
        samples = pd.DataFrame({'parent_name': ['parent1'],
                                'sample_name': ['sample1'],
                                'reference_file': ['ref1.fa']})
        # Should use all defaults if column is missing
        assert minimap2_params_with_default(wildcards, samples) == '-x map-hifi -B 1.5 --end-bonus 5'

    @pytest.mark.parametrize("invalid_val", [None, np.nan, 123, [], {}])
    def test_invalid_minimap2_params_value(self, invalid_val):
        wildcards = SimpleNamespace(sample='sample1')
        samples = pd.DataFrame({'parent_name': ['parent1'],
                                'sample_name': ['sample1'],
                                'reference_file': ['ref1.fa'],
                                'minimap2_params': [invalid_val]})
        # Should error if not string
        with pytest.raises(ValueError):
            minimap2_params_with_default(wildcards, samples)

    def test_sample_not_found(self):
        wildcards = SimpleNamespace(sample='not_in_samples')
        samples = pd.DataFrame({'parent_name': ['parent1'],
                                'sample_name': ['sample1'],
                                'reference_file': ['ref1.fa'],
                                'minimap2_params': ['']})
        with pytest.raises(ValueError):
            minimap2_params_with_default(wildcards, samples)

    def test_sample_is_both_sample_and_parent(self):
        wildcards = SimpleNamespace(sample='sample1')
        samples = pd.DataFrame({'parent_name': ['sample1'],
                                'sample_name': ['sample1'],
                                'reference_file': ['ref1.fa'],
                                'minimap2_params': ['']})
        with pytest.raises(ValueError):
            minimap2_params_with_default(wildcards, samples)

class TestGetParents:
    
    @pytest.mark.parametrize("parent, exp", [
        ('parent1', 'file1.fa'),
        ('parent2', 'file2.fa'),
        ('parent3', 'error'),
        ('not a parent', 'error'),
    ])
    def test_get_parents(self, parent, exp):

        wildcards = SimpleNamespace(sample=parent)
        samples = pd.DataFrame({'parent_name': ['parent1', 'parent2', 'parent3'], 
                                'parent_file': ['file1.fa', 'file2.fa', np.nan]})
        if exp == 'error':
            with pytest.raises(ValueError):
                get_parents(wildcards, samples)
        else:
            result = get_parents(wildcards, samples)
            assert result == exp

class TestFillParents:

    @pytest.mark.parametrize("sample, expected", [
        ('sample1', 'parent1.fa'),
        ('sample2', 'parent1.fa'),
        ('sample3', 'parent3.fa'),
        ('sample4', 'error'),
    ])
    def test_fill_parents(self, sample, expected):
        
        # set up
        wildcards = SimpleNamespace(sample=sample)
        filename = '{sample}.fa'
        samples = pd.DataFrame({'parent_name': ['parent1', 'parent1', 'parent3'], 
                                'sample_name': ['sample1', 'sample2', 'sample3']})
        
        # run fucntion and check
        if expected == 'error':
            with pytest.raises(KeyError):
                fill_parents(wildcards, samples, filename)
        else:
            result = fill_parents(wildcards, samples, filename)
            assert result == [expected]

class TestGetDmatInput:

    @pytest.mark.parametrize("seq_type, nt_seq, aa_seq, exp", [
        ("nt-seq", "nt.fa", "aa.fa", "nt.fa"),
        ("aa-seq", "nt.fa", "aa.fa", "aa.fa"),
        ('blah', "nt.fa", "aa.fa", 'error')
    ])
    def test_get_dmat_input_nt_seq(self, seq_type, nt_seq, aa_seq, exp):

        wildcards = SimpleNamespace(seq_type=seq_type)
        if exp == 'error':
            with pytest.raises(ValueError):
                get_dmat_input(wildcards, nt_seq, aa_seq)
        else:
            result = get_dmat_input(wildcards, "nt.fa", "aa.fa")
            assert result == exp

class TestGetReadsForCounting:

    @pytest.mark.parametrize("tech", ('np-cc', 'np', 'pb'))
    @pytest.mark.parametrize("read_file", ("file.fastq", "file.fastq.gz", "file.fq", "file.fq.gz"))
    @pytest.mark.parametrize("min_reps", (3, 1, None, np.nan))
    def test_get_reads_for_counting(self, tech, read_file, min_reps):

        # set up
        wildcards = SimpleNamespace(sample = "sample1")
        samples = pd.DataFrame({'sample_name': ["sample1"], 'seq_tech': [tech], 'read_file': [read_file], 'min_reps': [min_reps]})

        # run function
        res = get_reads_for_counting(wildcards, samples, "consensus.fa", "consesus_filt.fa")

        # check
        if tech != 'np-cc':
            assert res == [read_file]
        else:
            exp = [read_file, "consensus.fa"]
            if min_reps is not None and not np.isnan(min_reps):
                exp.append("consesus_filt.fa")
            assert res == exp

class TestFormatInputReads:

    @pytest.mark.parametrize("input, exp", [
        (['file.fastq',], '--fastq-files file.fastq'),
        (['file.fq',], '--fastq-files file.fq'),
        (['file.fa',], '--fasta-files file.fa '),
        (['file.fastq', 'file.fq'], '--fastq-files file.fastq file.fq'),
        (['file.fq', 'file.fa'], '--fasta-files file.fa --fastq-files file.fq'),
        (['file.fastq', 'file.fq', 'file.fa'], '--fasta-files file.fa --fastq-files file.fastq file.fq'),
    ])

    def test_format_input_reads(self, input, exp):
        
        assert format_input_reads(input) == exp
        