import pytest

from scripts.variant_frequency_long import get_variant_frequency

class TestGetVariantFrequency:

    @pytest.mark.parametrize("write_variants", 
                             [(True, 'some_variants'), 
                              (False, 'some_variants')], 
                            indirect=['write_variants'])
    def test_get_variant_frequency(self, write_variants):
        
        # get variants
        _, expected_vars, temp = write_variants
        expected_counts = {var: 1 for var in expected_vars}

        # get frequency
        counts = get_variant_frequency(temp.name)

        # check
        assert counts == expected_counts

        # clean up
        temp.close()

    @pytest.mark.parametrize("write_variants_repeated", 
                             [(True, (1,2,3,4,5,6), 'some_variants'), 
                              (False, (1,2,3,4,5,6), 'some_variants'),
                              (True, (2,2,2,2,2,2), 'some_variants'), 
                              (False, (2,2,2,2,2,2), 'some_variants')], 
                            indirect=['write_variants_repeated'])
    def test_get_variant_frequency(self, write_variants_repeated):
        
        # get variants
        _, n_repeats, expected_vars, temp = write_variants_repeated
        read_count = max(n_repeats)
        expected_counts = {var: rep/read_count for rep, var in zip(n_repeats, expected_vars)}

        # get frequency
        counts = get_variant_frequency(temp.name)

        # check
        assert counts == expected_counts

        # clean up
        temp.close()