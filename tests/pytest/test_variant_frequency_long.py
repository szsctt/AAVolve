import tempfile
import pytest

from scripts.variant_frequency_long import get_variant_frequency, write_freqs, write_variants, main

class TestGetVariantFrequency:

    @pytest.mark.parametrize("write_vars", 
                             [(True, 'some_variants'), 
                              (False, 'some_variants')], 
                            indirect=['write_vars'])
    def test_get_variant_frequency(self, write_vars):
        
        # get variants
        _, expected_vars, temp = write_vars
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


class TestWriteFreqs:

    @pytest.mark.parametrize("include_parents", [True, False]) 
    def test_write_freqs(self, some_variants, some_variant_frequencies, include_parents):
        
        if include_parents:
            parents = set(some_variants)
            par = "parental"
        else:
            parents = set()
            par = "non_parental"

        # write file
        with tempfile.NamedTemporaryFile(mode='w+') as temp:
            write_freqs(some_variant_frequencies, parents, temp.name)

            # read file
            temp.seek(0)
            lines = temp.readlines()

        # check
        assert len(lines) == len(some_variants) + 1
        assert lines[0] == 'query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n'
        for line, var in zip(lines[1:], some_variants):
            assert line == var.print_line(par)[:-1] + f"\t{some_variant_frequencies[var]}\n"
          

    def test_write_freqs_2(self, some_variants, some_variant_frequencies):
        
        parents = set(some_variants[:2])

        # write file
        with tempfile.NamedTemporaryFile(mode='w+') as temp:
            write_freqs(some_variant_frequencies, parents, temp.name)

            # read file
            temp.seek(0)
            lines = temp.readlines()

        # check
        assert len(lines) == len(some_variants) + 1
        assert lines[0] == 'query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n'
        for line, var in zip(lines[1:], some_variants):
            par = 'parental' if var in parents else 'non_parental'  
            assert line == var.print_line(par)[:-1] + f"\t{some_variant_frequencies[var]}\n"           
            
class TestWriteVariants:

    @pytest.mark.parametrize("write_vars", 
                             [(True, 'some_variants'), 
                              (False, 'some_variants')], 
                            indirect=['write_vars'])
    def test_write_variants(self, write_vars):
        
        shorter, expected_vars, temp = write_vars
        var_freqs = {var: 1 for var in expected_vars}

        # write file
        with tempfile.NamedTemporaryFile(mode='w+') as temp2:
            write_variants(var_freqs, temp2.name, temp.name)

            # read file
            temp2.seek(0)
            lines = temp2.readlines()

        # check
        assert len(lines) == len(expected_vars) + 1
        if shorter:
          assert lines[0] == 'query_name\tpos\tref_bases\tquery_bases\taa_change\n'
          ref = None
        else:
          assert lines[0] == 'reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n'
          ref = 'reference'
        for line, var in zip(lines[1:], expected_vars):
            assert line == var.print_line(query_name='non_parental', ref_name=ref)

        # clean up
        temp.close()
      

class TestMain:

    def test_main(self, resultfile_aav2389_some, resultfile_aav2389, monkeypatch):
        
        expected_high = ['reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n', 
                         'AAV2\t55\tnon_parental\tA56C\tA\tC\tFalse\n']
        expected_all = ['query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n', 
                        'parental\t40\tC\tA\tTrue\t1.0\n', 
                        'parental\t41\tT\tC\tTrue\t1.0\n', 
                        'parental\t44\tC\tT\tFalse\t1.0\n', 
                        'parental\t53\tA\tC\tFalse\t1.0\n', 
                        'non_parental\t55\tA\tC\tFalse\t0.6666666666666666\n']
        expected_par = ['query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n', 
                        'parental\t40\tC\tA\tTrue\t1.0\n', 
                        'parental\t41\tT\tC\tTrue\t1.0\n', 
                        'parental\t44\tC\tT\tFalse\t1.0\n', 
                        'parental\t53\tA\tC\tFalse\t1.0\n']
        
        with (tempfile.NamedTemporaryFile(mode='w+') as high,
              tempfile.NamedTemporaryFile(mode='w+') as all,
              tempfile.NamedTemporaryFile(mode='w+') as par):
            # set sys.argv
            monkeypatch.setattr('sys.argv', 
                                ['variant_frequency_long.py', 
                                 '-i', resultfile_aav2389_some,
                                 '-p', resultfile_aav2389, 
                                 '-o', high.name,
                                 '-oa', all.name,
                                  '-op', par.name,
                                  '-f', '0.5'])

            main()

            # read files
            high.seek(0), all.seek(0), par.seek(0)
            high_lines = high.readlines()
            all_lines = all.readlines()
            par_lines = par.readlines()
        
        assert high_lines == expected_high
        assert all_lines == expected_all
        assert par_lines == expected_par