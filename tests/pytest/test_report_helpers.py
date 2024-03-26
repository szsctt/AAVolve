import tempfile
import pytest
import numpy as np
import pandas as pd
import plotly.express as px

from scripts.report_helpers import (
    import_read_count_data, assign_file_type, read_count_graph, read_fraction_graph, print_fraction_nt_reads_pass,
    print_unique_nt_reads, print_unique_aa_reads, read_assigned_parents, parent_heatmap, plot_breakpoints,
    plot_parent_frequencies, make_distance_heatmap, parent_colors, numeric_position
    )

@pytest.fixture
def example_dfs():
    df_in = pd.DataFrame({'filename': ['file1.fastq', 'read_ids.txt', 'pivoted.tsv', 'file2-nt-seq-counts.tsv', 'file3-aa-seq-counts.tsv'], 
                          'file_type': ['fastq', 'variant_tsv', 'pivoted_tsv', 'distinct_read_counts', 'distinct_read_counts'], 
                          'Count': [200, 100, 99, 50, 10], 
                          })
    df_exp = pd.DataFrame({'filename': ['file1.fastq', 'read_ids.txt', 'pivoted.tsv', 'file2-nt-seq-counts.tsv', 'file3-aa-seq-counts.tsv'], 
                           'file_type': ['fastq', 'variant_tsv', 'pivoted_tsv', 'distinct_read_counts', 'distinct_read_counts'], 
                           'Count': [200, 100, 99, 50, 10], 
                           'File type': ['Input', 'Filtered by reference coverage', 'Filtered non-parental variants', 'Distinct at nucleotide level', 'Distinct at amino acid level'],
                           'Fraction of reads': [1.0, 0.5, 0.495, 0.25, 0.05]
                           })
    return df_in, df_exp

@pytest.fixture()
def example_parent_df():

    df = pd.DataFrame({'count': [100, 100, 100], '40:sub': ['AAV2', 'AAV2', 'AAV2,AAV3'], '41_42:del': ['AAV2', 'AAV2', 'AAV3']})

    return df

@pytest.fixture
def example_freq_df():

    df = pd.DataFrame({'variant': ['40:sub', '40:sub', '41_42:del', '41_42:del'], 'parent': ['AAV2', 'AAV3', 'AAV2', 'AAV3'], 'frequency': [0.5, 0.5, 0.5, 0.5]})

    return df

@pytest.fixture
def example_freq_df_longer():

    df = pd.DataFrame({'variant': ['40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub', '40:sub',
                                   '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del', '41_42:del',], 
                       'parent': ['AAV2', 'AAV3', 'AAV4', 'AAV5', 'AAV6', 'AAV7', 'AAV8', 'AAV9', 'AAV10', 'AAV12', 'LK03', 'SYD12', 'HRP5', 'rh74',
                                  'AAV2', 'AAV3', 'AAV4', 'AAV5', 'AAV6', 'AAV7', 'AAV8', 'AAV9', 'AAV10', 'AAV12', 'LK03', 'SYD12', 'HRP5', 'rh74',], 
                       'frequency': [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, ]})

    return df

@pytest.fixture
def example_freq_def_longest():

    df = pd.DataFrame({'variant': ['40:sub']*25,
                       'parent': [f'par{i}' for i in range(25)],
                       'frequency': [0.5]*25})
    
    return df


@pytest.fixture
def example_breakpoint_df():
    
    df = pd.DataFrame({'location': ['40:sub', '41_42:del'], 'breakpoints': [1, 0]})

    return df

@pytest.fixture
def example_dmat():


    dmat = np.array([[0, 0.1], [0.1, 0]])

    return dmat

class TestImportReadCountData:
    
    def test_import_read_count_data(self, example_dfs):
        
        # create a temporary file
        df_in, df_exp = example_dfs
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            df = import_read_count_data(f.name, 'np-cc')
        
        # check the output
        assert all(df == df_exp)

    def test_import_read_count_data_error_unknown_file_type(self):
        
        # create a temporary file
        df_in = pd.DataFrame({'filename': ['file1.fastq', 'read_ids.txt'], 'file_type': ['fastq', 'variants'], 'Count': [200, 100]})
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            with pytest.raises(ValueError):
                df = import_read_count_data(f.name, 'np-cc')

    def test_import_read_count_data_error_multiple_inputs(self):
        
        # create a temporary file
        df_in = pd.DataFrame({'filename': ['file1.fastq', 'file2.fastq'], 'file_type': ['fastq', 'fastq'], 'Count': [200, 100]})
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            with pytest.raises(AssertionError):
                df = import_read_count_data(f.name, 'np-cc')

    def test_import_read_count_data_error_increasing_reads(self):
        
        # create a temporary file
        df_in = pd.DataFrame({'filename': ['file1.fastq', 'read_ids.txt'], 'file_type': ['fastq', 'variants'], 'Count': [100, 200]})
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            with pytest.raises(ValueError):
                df = import_read_count_data(f.name, 'np-cc')

class TestAssignFileType:

    @pytest.mark.parametrize("file_name, file_type, seq_type, exp", [ 
        ('file1.fastq', 'fastq', 'np-cc', 'Input'),
        ('out/c3poa/sample/split/file1.fasta', 'fasta', 'np-cc', 'Consensus'),
        ('out/c3poa_filt/file1.fasta', 'fasta', 'np-cc', 'Filtered by repeats'),
        ('file1.fasta', 'fasta', 'sg', 'Input'),
        ('file1.txt', 'variant_tsv', 'np-cc', 'Filtered by reference coverage'),
        ('file1.tsv', 'pivoted_tsv', 'np-cc', 'Filtered non-parental variants'),
        ('file1-parent-counts.tsv', 'distinct_read_counts', 'np-cc', 'Distinct at nucleotide level'),
        ('file1-nt-seq-counts.tsv', 'distinct_read_counts', 'np-cc', 'Distinct at nucleotide level'),
        ('file1-aa-seq-counts.tsv', 'distinct_read_counts', 'np-cc', 'Distinct at amino acid level'),
        ('file1.tsv', 'pivoted_tsv', 'np', 'Filtered non-parental variants'),
        ('file1-nt-seq-counts.tsv', 'distinct_read_counts', 'np', 'Distinct at nucleotide level'),
        ('file1-aa-seq-counts.tsv', 'distinct_read_counts', 'np', 'Distinct at amino acid level'),
    ])
    def test_assign_file_type(self, file_name, file_type, seq_type, exp):

        # test the function
        res = assign_file_type(file_name, file_type, seq_type)

        # check the output
        assert res == exp


    @pytest.mark.parametrize("file_name, file_type, seq_type", [ 
        ('file1.fasta', 'fasta', 'unknown'),
        ('file1.txt', 'unknown', 'np-cc'),
        ('file1.txt', 'unknown', 'unknown'),
    ])
    def test_assign_file_type_error_unknown_file_type(self, file_name, file_type, seq_type):

        # test the function
        with pytest.raises(ValueError):
            res = assign_file_type(file_name, file_type, seq_type)

class TestReadCountGraph:

    def test_read_count_graph(self, example_dfs):

        # create a temporary file
        df_in, df_exp = example_dfs
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            fig = read_count_graph(f.name, 'np-cc')

        # generate expected figure
        df_exp = df_exp[df_exp['File type'] != "Distinct at nucleotide level"]
        df_exp = df_exp[df_exp['File type'] != "Distinct at amino acid level"]
        fig_exp = px.line(df_exp, x='File type', y='Count', markers=True,
                  color_discrete_sequence=['black'])
        fig_exp.update_xaxes(tickangle=90)

        # check the output
        assert fig == fig_exp

class TestReadFractionGraph:

    def test_read_fraction_graph(self, example_dfs):

        # create a temporary file
        df_in, df_exp = example_dfs
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            fig = read_fraction_graph(f.name, 'np-cc')

        # generate expected figure
        df_exp = df_exp[df_exp['File type'] != "Distinct at nucleotide level"]
        df_exp = df_exp[df_exp['File type'] != "Distinct at amino acid level"]
        fig_exp = px.line(df_exp, x='File type', y='Fraction of reads', markers=True, 
                  color_discrete_sequence=['black'])
        fig_exp.update_xaxes(tickangle=90)

        # check the output
        assert fig == fig_exp


class TestPrintFractionNtReadsPass:
    def test_print_fraction_nt_reads_pass(self, example_dfs, capsys):

        # create a temporary file
        df_in, df_exp = example_dfs
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            print_fraction_nt_reads_pass(f.name, 'np-cc')

        # check the output
        pass_frac = df_exp[df_exp['File type'] == "Distinct at nucleotide level"]['Fraction of reads'].to_list()[0]
        captured = capsys.readouterr()
        assert captured.out == f'{pass_frac*100:.2f}%\n'

class TestPrintUniqueNtReads:

    def test_print_unique_nt_reads(self, example_dfs, capsys):

        # create a temporary file
        df_in, df_exp = example_dfs
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            print_unique_nt_reads(f.name, 'np-cc')

        # check the output
        count = df_exp[df_exp['File type'] == "Distinct at nucleotide level"]['Count'].to_list()[0]
        captured = capsys.readouterr()
        assert captured.out == f'{count}\n'


class TestPrintUniqueAaReads:

    def test_print_unique_aa_reads(self, example_dfs, capsys):

        # create a temporary file
        df_in, df_exp = example_dfs
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False, header=False)
            f.seek(0)

            # test the function
            print_unique_aa_reads(f.name, 'np-cc')

        # check the output
        count = df_exp[df_exp['File type'] == "Distinct at amino acid level"]['Count'].to_list()[0]
        captured = capsys.readouterr()
        assert captured.out == f'{count}\n'

class TestReadAssignedParents:

    def test_read_assigned_parents(self, example_parent_df):

        # create a temporary file
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            example_parent_df.to_csv(f.name, sep='\t', index=False)
            f.seek(0)

            # test the function
            df = read_assigned_parents(f.name)

        # check the output
        assert all(df == example_parent_df)

    
    def test_read_assigned_parents_max(self, example_parent_df):

        # make a dataframe with more than 1000 rows
        df_in = pd.concat([example_parent_df]*1000, ignore_index=True)

        # create a temporary file
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            df_in.to_csv(f.name, sep='\t', index=False)
            f.seek(0)

            # test the function
            df = read_assigned_parents(f.name)

        # check the output
        assert all(df == df_in.iloc[:1000]) 

class TestParentHeatmap:

    def test_parent_heatmap(self, example_parent_df, example_freq_df):

        exp_fig = "Figure({\n    'data': [{'line': {'color': 'black'},\n              'marker': {'color': 'black'},\n              'mode': 'lines+markers',\n              'showlegend': False,\n              'type': 'scatter',\n              'x': array([100, 100, 100]),\n              'xaxis': 'x',\n              'y': array([0, 1, 2]),\n              'yaxis': 'y'},\n             {'colorbar': {'ticktext': [AAV2, AAV3, non parental, multiple],\n                           'tickvals': [0.125, 0.375, 0.625, 0.875]},\n              'colorscale': [[0.0, '#636EFA'], [0.25, '#636EFA'], [0.25,\n                             '#EF553B'], [0.5, '#EF553B'], [0.5, '#00CC96'], [0.75,\n                             '#00CC96'], [0.75, '#AB63FA'], [1.0, '#AB63FA']],\n              'hovertemplate': '{text}<extra></extra>',\n              'text': array([['AAV2', 'AAV2'],\n                             ['AAV2', 'AAV2'],\n                             ['multiple', 'AAV3']], dtype=object),\n              'type': 'heatmap',\n              'x': array(['40', '41_42'], dtype=object),\n              'xaxis': 'x2',\n              'yaxis': 'y2',\n              'z': array([[0.  , 0.  ],\n                          [0.  , 0.  ],\n                          [0.75, 0.25]]),\n              'zmax': 1,\n              'zmin': 0}],\n    'layout': {'template': '...',\n               'xaxis': {'anchor': 'y', 'domain': [0.0, 0.19], 'title': {'text': 'Count'}},\n               'xaxis2': {'anchor': 'y2', 'domain': [0.24, 1.0], 'title': {'text': 'Position in reference'}},\n               'yaxis': {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': 'Read'}},\n               'yaxis2': {'anchor': 'x2', 'domain': [0.0, 1.0], 'matches': 'y', 'showticklabels': False}}\n})"

        # create a temporary file
        with (tempfile.NamedTemporaryFile(mode='w+') as f,
              tempfile.NamedTemporaryFile(mode='w+') as f2):
            example_parent_df.to_csv(f.name, sep='\t', index=False)
            example_freq_df.to_csv(f2.name, sep='\t', index=False)
            f.seek(0), f2.seek(0)

            # test the function
            fig = parent_heatmap(f.name, f2.name)

        # check the output
        assert fig.__str__() == exp_fig

class TestPlotBreakpoints:

    def test_plot_breakpoints(self, example_dfs, example_breakpoint_df):

        exp_fig = "Figure({\n    'data': [{'hovertemplate': 'Position in reference=%{x}<br>Breakpoint frequency (%)=%{y}<extra></extra>',\n              'legendgroup': '',\n              'line': {'color': 'black', 'dash': 'solid'},\n              'marker': {'symbol': 'circle'},\n              'mode': 'lines',\n              'name': '',\n              'orientation': 'v',\n              'showlegend': False,\n              'type': 'scatter',\n              'x': array(['40', '41_42  '], dtype=object),\n              'xaxis': 'x',\n              'y': array([1.01010101, 0.        ]),\n              'yaxis': 'y'}],\n    'layout': {'legend': {'tracegroupgap': 0},\n               'margin': {'t': 60},\n               'template': '...',\n               'xaxis': {'anchor': 'y', 'domain': [0.0, 1.0], 'title': {'text': 'Position in reference'}},\n               'yaxis': {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': 'Breakpoint frequency (%)'}}}\n})" 

        # create a temporary file
        df_in = example_dfs[0]

        with (tempfile.NamedTemporaryFile('w+') as breakf, 
              tempfile.NamedTemporaryFile('w+') as countf):
            df_in.to_csv(countf.name, sep='\t', index=False, header=False)
            example_breakpoint_df.to_csv(breakf.name, sep='\t', index=False)

            fig = plot_breakpoints(breakf.name, countf.name, 'np-cc')

        assert str(fig) == exp_fig

class TestPlotParentFrequencies:

    def test_plot_parent_frequencies(self, example_freq_df):

        exp_fig =  "Figure({\n    'data': [{'alignmentgroup': 'True',\n              'hovertemplate': ('Parent=AAV2<br>Position in ref' ... 'quency (%)=%{y}<extra></extra>'),\n              'legendgroup': 'AAV2',\n              'marker': {'color': '#636EFA', 'line': {'width': 0}, 'pattern': {'shape': ''}},\n              'name': 'AAV2',\n              'offsetgroup': 'AAV2',\n              'orientation': 'v',\n              'showlegend': True,\n              'textposition': 'auto',\n              'type': 'bar',\n              'x': array(['40', '41_42  '], dtype=object),\n              'xaxis': 'x',\n              'y': array([50., 50.]),\n              'yaxis': 'y'},\n             {'alignmentgroup': 'True',\n              'hovertemplate': ('Parent=AAV3<br>Position in ref' ... 'quency (%)=%{y}<extra></extra>'),\n              'legendgroup': 'AAV3',\n              'marker': {'color': '#EF553B', 'line': {'width': 0}, 'pattern': {'shape': ''}},\n              'name': 'AAV3',\n              'offsetgroup': 'AAV3',\n              'orientation': 'v',\n              'showlegend': True,\n              'textposition': 'auto',\n              'type': 'bar',\n              'x': array(['40', '41_42  '], dtype=object),\n              'xaxis': 'x',\n              'y': array([50., 50.]),\n              'yaxis': 'y'}],\n    'layout': {'bargap': 0,\n               'bargroupgap': 0,\n               'barmode': 'relative',\n               'legend': {'title': {'text': 'Parent'}, 'tracegroupgap': 0},\n               'margin': {'t': 60},\n               'template': '...',\n               'xaxis': {'anchor': 'y',\n                         'categoryarray': array(['40', '40', '41_42  ', '41_42  '], dtype=object),\n                         'categoryorder': 'array',\n                         'domain': [0.0, 1.0],\n                         'title': {'text': 'Position in reference'}},\n               'yaxis': {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': 'Parent frequency (%)'}}}\n})"

        with tempfile.NamedTemporaryFile('w+') as f:

            example_freq_df.to_csv(f.name, sep='\t', index=False)

            fig = plot_parent_frequencies(f.name)

        assert str(fig) == exp_fig

class TestMakeDistanceHeatmap:

    def test_make_distance_heatmap(self, example_dmat):

        exp_fig = "Figure({\n    'data': [{'type': 'heatmap', 'z': array([[0. , 0.1],\n          [0.1, 0. ]])}], 'layout': {'template': '...'}\n})"

        with tempfile.NamedTemporaryFile('w+') as f:

            np.savetxt(f.name, example_dmat)

            fig = make_distance_heatmap(f.name)

        assert str(fig) == exp_fig


class TestParentColors:

    def test_parent_colors(self, example_freq_df):

        exp_dict = {'AAV2': '#636EFA', 'AAV3': '#EF553B', 'non parental': '#00CC96', 'multiple': '#AB63FA'}

        with tempfile.NamedTemporaryFile('w+') as f:

            example_freq_df.to_csv(f.name, sep='\t', index=False)

            res = parent_colors(f.name)

        assert res == exp_dict

    def test_parent_colors_2(self, example_freq_df_longer):

        exp_dict = {'AAV2': '#FD3216', 'AAV3': '#00FE35', 'AAV4': '#6A76FC', 'AAV5': '#FED4C4', 'AAV6': '#FE00CE', 'AAV7': '#0DF9FF', 'AAV8': '#F6F926', 'AAV9': '#FF9616', 'AAV10': '#479B55', 'AAV12': '#EEA6FB', 'LK03': '#DC587D', 'SYD12': '#D626FF', 'HRP5': '#6E899C', 'rh74': '#00B5F7', 'non parental': '#B68E00', 'multiple': '#C9FBE5'}

        with tempfile.NamedTemporaryFile('w+') as f:

            example_freq_df_longer.to_csv(f.name, sep='\t', index=False)

            res = parent_colors(f.name)

        assert res == exp_dict


    def test_parent_colors_longest(self, example_freq_def_longest):

        exp_dict = {'par0': 'rgb(48, 18, 59)', 'par1': 'rgb(57, 45, 119)', 'par2': 'rgb(65, 73, 176)', 'par3': 'rgb(68, 99, 212)', 'par4': 'rgb(68, 124, 239)', 'par5': 'rgb(61, 148, 247)', 'par6': 'rgb(50, 172, 243)', 'par7': 'rgb(34, 197, 221)', 'par8': 'rgb(30, 216, 198)', 'par9': 'rgb(35, 232, 173)', 'par10': 'rgb(59, 242, 144)', 'par11': 'rgb(92, 251, 112)', 'par12': 'rgb(128, 252, 85)', 'par13': 'rgb(164, 252, 59)', 'par14': 'rgb(188, 241, 55)', 'par15': 'rgb(212, 229, 52)', 'par16': 'rgb(230, 211, 56)', 'par17': 'rgb(245, 191, 56)', 'par18': 'rgb(251, 168, 49)', 'par19': 'rgb(251, 142, 39)', 'par20': 'rgb(246, 112, 27)', 'par21': 'rgb(235, 86, 16)', 'par22': 'rgb(221, 63, 8)', 'par23': 'rgb(202, 44, 4)', 'par24': 'rgb(180, 27, 1)', 'non parental': 'rgb(152, 15, 1)', 'multiple': 'rgb(122, 4, 2)'}

        with tempfile.NamedTemporaryFile('w+') as f:

            example_freq_def_longest.to_csv(f.name, sep='\t', index=False)

            res = parent_colors(f.name)

        assert res == exp_dict


class TestNumericPosition:

    @pytest.mark.parametrize("variant, exp", [
            ("40:sub", "40"),
            ("40_41:del", "40_41  "), 
            ("40:ins", "40 ")
    ])
    def test_numeric_position(self, variant, exp):

        res = numeric_position(pd.Series(variant))

        assert all(res == pd.Series(exp))