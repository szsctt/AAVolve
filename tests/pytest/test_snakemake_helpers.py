from types import SimpleNamespace
import pytest
import numpy as np
import pandas as pd
from scripts.snakemake_helpers import (
    get_column_by_sample, is_fastq, get_reads_for_counting, get_dmat_input, format_input_reads, get_parents, fill_parents, get_reference
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
        