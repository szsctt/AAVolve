import os
import plotly_express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from utils import MAX_SEQS

#### counts of reads ####

def import_read_count_data(df_file, seq_type):

    df = pd.read_csv(df_file, delimiter='\t', header=None, names=['filename', 'file_type', 'Count'])

    df['Fraction of reads'] = df['Count'] / df['Count'][0]

    df['File type'] = df.apply(lambda row: assign_file_type(row['filename'], row['file_type'], seq_type), axis=1)
   
    return df

def assign_file_type(file_name, file_type, seq_type):
    if file_type == "fastq":
        if seq_type == "np-cc":
            dir = os.path.dirname(file_name)
            dir2 = os.path.dirname(dir)
            dir3 = os.path.dirname(dir2)
            if os.path.basename(dir3) == "c3poa":
                return "Consensus"
            elif os.path.basename(dir) == "c3poa_filt":
                return "Filtered by repeats"
        else:
            return 'Input'
    elif file_type == "variant_tsv":
        return "Filtered by reference coverage"
    elif file_type == "pivoted_tsv":
        return "Filtered non-parental variants"
    elif file_type == "distinct_read_counts":
        if 'parent-counts' in os.path.basename(file_name):
            return "Distinct at nucleotide level"
        elif 'nt-seq-counts' in os.path.basename(file_name):
            return "Distinct at nucleotide level"
        elif 'aa-seq-counts' in os.path.basename(file_name):
            return "Distinct at amino acid level"

def read_count_graph(df_file, seq_type):

    df = import_read_count_data(df_file, seq_type)
    df = df[df['File type'] != "Distinct at nucleotide level"]
    df = df[df['File type'] != "Distinct at amino acid level"]
            

    fig = px.line(df, x='File type', y='Count', markers=True,
                  color_discrete_sequence=['black'])
    fig.update_xaxes(tickangle=90)

    return fig

def read_fraction_graph(df_file, seq_type):

    df = import_read_count_data(df_file, seq_type)
    df = df[df['File type'] != "Distinct at nucleotide level"]
    df = df[df['File type'] != "Distinct at amino acid level"]

    fig = px.line(df, x='File type', y='Fraction of reads', markers=True, 
                  color_discrete_sequence=['black'])
    fig.update_xaxes(tickangle=90)

    return fig

def print_fraction_nt_reads_pass(df_file, seq_type):

    df = import_read_count_data(df_file, seq_type)

    pass_frac = df[df['File type'] == "Distinct at nucleotide level"]['Fraction of reads'].to_list()[0]

    print(f'{pass_frac*100:.2f}%')

def print_unique_nt_reads(df_file, seq_type):

    df = import_read_count_data(df_file, seq_type)

    # get value for count
    count = df[df['File type'] == "Distinct at nucleotide level"]['Count'].to_list()[0]

    print(count)

def print_unique_aa_reads(df_file, seq_type):

    df = import_read_count_data(df_file, seq_type)

    # get value for count
    count = df[df['File type'] == "Distinct at amino acid level"]['Count'].to_list()[0]

    print(count)


#### variants ####
    

def read_assigned_parents(filename):

    df = pd.read_csv(filename, delimiter='\t', nrows=MAX_SEQS)

    return df


def parent_heatmap(filename, parent_freq_file):

    # https://chart-studio.plotly.com/~empet/15229/heatmap-with-a-discrete-colorscale/#/

    # read in data
    df = read_assigned_parents(filename)

    # reverse order of rows to get most frequent at top
    df = df.iloc[::-1]

    # get counts
    counts = df['count']
    
    # get a list of unique parents
    df = (
        df
          # drop count column
          .drop(columns=['count']) 
          # replace non_parental_1 etc with non_parental
          .apply(lambda x: x.str.replace("non_parental_\d+", "non parental", regex=True)) 
          # replace any values with commas with "multiple"
          .apply(lambda x: [i if ',' not in i else 'multiple' for i in x])
          ) 

    # get colors for each parent
    color_dict = parent_colors(parent_freq_file)
    parents = list(color_dict.keys())
    colors = list(color_dict.values())

    # convert df to numeric, with numbers between 
    # 0 and 1 corresponding to parents
    lnsp = np.linspace(0, 1, len(parents)+1)
    conv = dict(zip(parents, lnsp))
    df_nums = df.applymap(lambda x: conv[x])
    
    # determines which values are mapped to each color
    colorsc = []
    for i in range(len(colors)):
        colorsc.append([lnsp[i], colors[i]])
        colorsc.append([lnsp[i+1], colors[i]])

    
    # tick locations should be between changes in colorbar
    tickvals = [np.mean(lnsp[i:i+2]) for i in range(len(lnsp)-1)]

    # x axis should just be numbers, remove :sub etc
    df_nums.columns = df_nums.columns.str.replace(":.*", "", regex=True)
    p1 = go.Heatmap(
            z=df_nums.values, 
            x=numeric_position(df_nums.columns), 
            colorscale=colorsc, 
            zmin=0, zmax=1,
            hovertemplate='{text}<extra></extra>',
            text = df.values,
            colorbar=dict(tickvals=tickvals, ticktext=parents))

    p2 = go.Scatter(x=counts, 
                    y=df_nums.index, 
                    mode='lines+markers', 
                    line_color='black', 
                    marker=dict(color='black'), 
                    showlegend=False)
    fig = make_subplots(rows=1, cols=2, column_widths=[0.2, 0.8], shared_yaxes=True, horizontal_spacing=0.05, vertical_spacing=0.05)
    fig.add_trace(p2, row=1, col=1)
    fig.add_trace(p1, row=1, col=2)
    fig['layout']['xaxis']['title'] = 'Count'
    fig['layout']['yaxis']['title'] = 'Read'
    fig['layout']['xaxis2']['title'] = 'Position in reference'


    return fig


def plot_breakpoints(breakpoints_file, counts_file):

    # get total reads passing all filters
    counts = import_read_count_data(counts_file, 'np')
    total_nt_reads = counts[counts['File type'] == 'Filtered non-parental variants']['Count'].to_list()[0]

    # import breakpoint counts by position
    df = pd.read_csv(breakpoints_file, delimiter='\t')
    
    # Normalize to total number of reads
    df['breakpoints'] = df['breakpoints'] / total_nt_reads * 100
    df['location'] = numeric_position(df['location'])

    # make plot
    fig = px.line(df, x='location', y='breakpoints',
                  labels = {'breakpoints': 'Breakpoint frequency (%)',
                            'location': 'Position in reference'},
                            color_discrete_sequence=['black']
    )
    return fig

def plot_parent_frequencies(parents_file):

    # read data
    df = pd.read_csv(parents_file, delimiter='\t')

    # change 'non_parental_1' etc to 'non parental'
    df['parent'] = df['parent'].astype(str).str.replace("non_parental_\d+", "non parental", regex=True)

    # convert frequency to percentage
    df['frequency'] = df['frequency'] * 100

    # TODO: deal with different kinds of variants at same location
    df['variant'] = numeric_position(df['variant'])

    # map parent names to colors
    parent_colors_dict = parent_colors(parents_file)

    fig = px.bar(df, x='variant', y='frequency', color='parent',
                 labels = {'frequency': 'Parent frequency (%)',
                           'variant': 'Position in reference',
                           'parent': 'Parent'},
                 color_discrete_map=parent_colors_dict,
                           )
    
    # remove white lines around bars
    # https://stackoverflow.com/questions/69553283/how-to-remove-white-lines-around-bars-using-plotly-express-bar
    fig.update_traces(marker_line_width = 0,
                  selector=dict(type="bar"))

    fig.update_layout(bargap=0,
                  bargroupgap = 0,
                 )

    return fig

def make_distance_heatmap(distance_file):

    dmat = np.loadtxt(distance_file)

    p = go.Heatmap(z=dmat)
    fig = go.Figure(data=p)
    return fig

def parent_colors(parents_file):

    # read data
    df = pd.read_csv(parents_file, delimiter='\t')

    # change 'non_parental_1' etc to 'non parental'
    df['parent'] = df['parent'].astype(str).str.replace("non_parental_\d+", "non parental", regex=True)

    # get list of unique parents
    parents = df['parent'].unique().tolist()

    # move non_parental and multiple to the end
    if 'non parental' in parents:
        parents.remove('non parental')
    parents.append('non parental')
    if 'multiple' in parents:
        parents.remove('multiple')
    parents.append('multiple')

    # create colormap for heatmap
    if len(parents) < 11:
        # for 10 or fewer parents, use safe qualitative colors
        colors = px.colors.qualitative.Plotly[:len(parents)]
    elif len(parents) < 25:
        # for 11-24 parents, use light24 colors
        colors = px.colors.qualitative.Light24[:len(parents)]
    else:
        # sample colors from a continuous colormap
        # this probably isn't going to look good
        colors = px.colors.sample_colorscale(
            px.colors.sequential.Turbo, len(parents))
        
    return dict(zip(parents, colors))



def numeric_position(col):

    col = col.astype(str).str.replace(":sub", "", regex=True)
    col = col.str.replace(":ins", " ", regex=True)
    col = col.str.replace(":del", "  ", regex=True)

    return col
