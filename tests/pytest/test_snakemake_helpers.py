from types import SimpleNamespace
import pytest
import pandas as pd
from scripts.snakemake_helpers import get_column_by_sample


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