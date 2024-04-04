import tempfile
import pytest

from aavolve.variant_frequency_long import get_variant_frequency, write_freqs, write_variants, main


class TestGetVariantFrequency:

    @pytest.mark.parametrize("write_vars", 
                             [(True, 'some_variants'), 
                              (False, 'some_variants')], 
                            indirect=['write_vars'])
    def test_get_variant_frequency(self, write_vars):
        
        # get variants
        _, expected_vars, temp_vars, temp_rids = write_vars
        expected_counts = {var: 1 for var in expected_vars}

        # get frequency
        counts = get_variant_frequency(temp_vars.name, temp_rids.name)

        # check
        assert counts == expected_counts

    @pytest.mark.parametrize("write_variants_repeated,extra_reads", 
                             [
                              ((True, (1,2,3,4,5,6), 'some_variants'), 0), 
                              ((False, (1,2,3,4,5,6), 'some_variants'), 0),
                              ((True, (2,2,2,2,2,2), 'some_variants'), 0), 
                              ((False, (2,2,2,2,2,2), 'some_variants'), 0),
                              ((True, (1,2,3,4,5,6), 'some_variants'), 2), 
                              ((False, (1,2,3,4,5,6), 'some_variants'), 2),
                              ((True, (2,2,2,2,2,2), 'some_variants'), 2), 
                              ((False, (2,2,2,2,2,2), 'some_variants'), 2),
                              ], 
                            indirect=['write_variants_repeated'])
    def test_get_variant_frequency_2(self, write_variants_repeated, extra_reads):
        
        # get variants
        _, n_repeats, expected_vars, temp_vars, temp_rids = write_variants_repeated
        read_count = max(n_repeats) + extra_reads
        expected_counts = {var: rep/read_count for rep, var in zip(n_repeats, expected_vars)}

        # write extra reads to file
        with open(temp_rids.name, 'a') as f:
            for i in range(extra_reads):
                f.write(f"extra_{i}\n")
        temp_rids.seek(0)

        # get frequency
        counts = get_variant_frequency(temp_vars.name, temp_rids.name)

        # check
        assert counts == expected_counts


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
        
        shorter, expected_vars, temp_vars, _ = write_vars
        var_freqs = {var: 1 for var in expected_vars}

        # write file
        with tempfile.NamedTemporaryFile(mode='w+') as temp2:
            write_variants(var_freqs, temp2.name, temp_vars.name)

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
        i = 0
        for line, var in zip(lines[1:], expected_vars):
            assert line == var.print_line(query_name=f'non_parental_{i}', ref_name=ref)
            i += 1
      
    def test_no_variants(self):
        
        # write file
        with tempfile.NamedTemporaryFile(mode='w+') as temp, tempfile.NamedTemporaryFile(mode='w+') as temp2:
            # write header to temp2
            temp2.write('query_name\tpos\tref_bases\tquery_bases\taa_change\n')
            temp2.seek(0)
            write_variants({}, temp.name, temp2.name)

            # read file
            temp.seek(0)
            lines = temp.readlines()

        # check
        assert len(lines) == 1
        assert lines[0] == 'query_name\tpos\tref_bases\tquery_bases\taa_change\n'


class TestMain:

    def test_main(self, resultfile_aav2389_some, resultfile_aav2389, monkeypatch):
        
        expected_high = ['reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n', 
                         'AAV2\t55\tnon_parental_0\tA56C\tA\tC\tFalse\n']
        expected_all = ['query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n', 
                        'non_parental\t40\tC\tG\tTrue\t0.3333333333333333\n', 
                        'parental\t41\tT\tC\tTrue\t1.0\n', 
                        'parental\t44\tC\tT\tFalse\t1.0\n', 
                        'parental\t53\tA\tC\tFalse\t1.0\n', 
                        'non_parental\t55\tA\tC\tFalse\t0.6666666666666666\n',
                        'parental\t40\tC\tA\tTrue\t0.6666666666666666\n']
        expected_par = ['query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n', 
                        'parental\t41\tT\tC\tTrue\t1.0\n', 
                        'parental\t44\tC\tT\tFalse\t1.0\n', 
                        'parental\t53\tA\tC\tFalse\t1.0\n',
                        'parental\t40\tC\tA\tTrue\t0.6666666666666666\n']
        
        with (tempfile.NamedTemporaryFile(mode='w+') as high,
              tempfile.NamedTemporaryFile(mode='w+') as all,
              tempfile.NamedTemporaryFile(mode='w+') as par):
            # set sys.argv
            monkeypatch.setattr('sys.argv', 
                                ['variant_frequency_long.py', 
                                '-i', resultfile_aav2389_some[0],
                                '-r', resultfile_aav2389_some[1],
                                 '-p', resultfile_aav2389[0], 
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

    def test_main_headers_only(self, monkeypatch):
            
            with (
                tempfile.NamedTemporaryFile(mode='w+') as reads_in,
                tempfile.NamedTemporaryFile(mode='w+') as read_ids_in,
                tempfile.NamedTemporaryFile(mode='w+') as par_in,
                tempfile.NamedTemporaryFile(mode='w+') as high_out,
                tempfile.NamedTemporaryFile(mode='w+') as all_out,
                tempfile.NamedTemporaryFile(mode='w+') as par_out):

                # write headers
                for f in (reads_in, par_in):
                    f.write('query_name\tpos\tref_bases\tquery_bases\taa_change\n')
                    f.seek(0)

                # set sys.argv
                monkeypatch.setattr('sys.argv', 
                ['variant_frequency_long.py', 
                    '-i', reads_in.name,
                    '-r', read_ids_in.name,
                    '-p', par_in.name, 
                    '-o', high_out.name,
                    '-oa', all_out.name,
                    '-op', par_out.name,
                    '-f', '0.5'])
                
                # run function
                main()
    
                # read files
                high_out.seek(0), all_out.seek(0), par_out.seek(0)
                high_lines = high_out.readlines()
                all_lines = all_out.readlines()
                par_lines = par_out.readlines()
            
            assert high_lines == ['query_name\tpos\tref_bases\tquery_bases\taa_change\n']
            assert all_lines == []
            assert par_lines == []

    def test_main_no_parents(self, resultfile_aav2389_some, monkeypatch):
            
        expected_high = [
            'reference_name\tpos\tquery_name\tvar\tref_bases\tquery_bases\taa_change\n', 
            'AAV2\t41\tnon_parental_0\tT42C\tT\tC\tTrue\n', 
            'AAV2\t44\tnon_parental_1\tC45T\tC\tT\tFalse\n', 
            'AAV2\t53\tnon_parental_2\tA54C\tA\tC\tFalse\n', 
            'AAV2\t55\tnon_parental_3\tA56C\tA\tC\tFalse\n',
            'AAV2\t40\tnon_parental_4\tC41A\tC\tA\tTrue\n', 
            ]
        expected_all = [
            'query_name\tpos\tref_bases\tquery_bases\taa_change\tfreq\n', 
            'non_parental\t40\tC\tG\tTrue\t0.3333333333333333\n', 
            'non_parental\t41\tT\tC\tTrue\t1.0\n', 
            'non_parental\t44\tC\tT\tFalse\t1.0\n', 
            'non_parental\t53\tA\tC\tFalse\t1.0\n', 
            'non_parental\t55\tA\tC\tFalse\t0.6666666666666666\n',
             'non_parental\t40\tC\tA\tTrue\t0.6666666666666666\n', 
            ]
        
        with (tempfile.NamedTemporaryFile(mode='w+') as high,
              tempfile.NamedTemporaryFile(mode='w+') as all):
            # set sys.argv
            monkeypatch.setattr('sys.argv', 
                                ['variant_frequency_long.py', 
                                '-i', resultfile_aav2389_some[0],
                                '-r', resultfile_aav2389_some[1],
                                '-o', high.name,
                                '-oa', all.name,
                                '-f', '0.5'])

            main()

            # read files
            high.seek(0), all.seek(0)
            high_lines = high.readlines()
            all_lines = all.readlines()
        
        assert high_lines == expected_high
        assert all_lines == expected_all