import os
import plotly_express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

MAX_ROWS = 10000

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
            

    fig = px.line(df, x='File type', y='Count', markers=True)
    fig.update_xaxes(tickangle=90)

    return fig

def read_fraction_graph(df_file, seq_type):

    df = import_read_count_data(df_file, seq_type)
    df = df[df['File type'] != "Distinct at nucleotide level"]
    df = df[df['File type'] != "Distinct at amino acid level"]

    fig = px.line(df, x='File type', y='Fraction of reads', markers=True)
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

    df = pd.read_csv(filename, delimiter='\t', nrows=MAX_ROWS)

    return df


def parent_heatmap(filename):

    # https://chart-studio.plotly.com/~empet/15229/heatmap-with-a-discrete-colorscale/#/

    # read in data
    df = read_assigned_parents(filename)
    counts = df['count']
    df = (
        df
          # drop count column
          .drop(columns=['count']) 
          # replace non_parental_1 etc with non_parental
          .apply(lambda x: x.str.replace("non_parental_\d+", "non_parental", regex=True)) 
          # replace any values with commans with "multiple"
          .apply(lambda x: [i if ',' not in i else 'multiple' for i in x])
          ) 
    df_long = df.melt()

    # get a list of unique parents
    parents = list(set(df_long.value.str.split(',').explode()))
    # move non_parental and multiple to the end
    if 'non_parental' in parents:
        parents.remove('non_parental')
        parents.append('non_parental')
    if 'multiple' in parents:
        parents.remove('multiple')
        parents.append('multiple')

    # convert df to numeric, with numbers between 
    # 0 and 1 corresponding to parents
    lnsp = np.linspace(0, 1, len(parents))
    conv = dict(zip(parents, lnsp))
    df_nums = df.applymap(lambda x: conv[x])

    # create colormap for heatmap
    if len(parents) < 11:
        # for 10 or fewer parents, use safe qualitative colors
        colors = px.colors.qualitative.Safe[:len(parents)]
    elif len(parents) < 25:
        # for 11-24 parents, use light24 colors
        colors = px.colors.qualitative.Light24[:len(parents)]
    else:
        # sample colors from a continuous colormap
        # this probably isn't going to look good
        colors = px.colors.sample_colorscale(
            px.colors.sequential.Turbo, lnsp)
    
    # determines which values are mapped to each color
    colorsc = [[i,j] for i,j in zip(lnsp, colors)]

    # x axis should just be numbers, remove :sub etc
    df_nums.columns = df_nums.columns.str.replace(":.*", "", regex=True)
    p1 = go.Heatmap(
            z=df_nums.values, 
            x=df_nums.columns, 
            colorscale=colorsc, 
            zmin=0, zmax=1,
            colorbar=dict(tickvals=lnsp, ticktext=parents))

    p2 = go.Scatter(x=counts, y=df_nums.index, mode='lines+markers')
    fig = make_subplots(rows=1, cols=2, column_widths=[0.2, 0.8], shared_yaxes=True)
    fig.add_trace(p2, row=1, col=1)
    fig.add_trace(p1, row=1, col=2)


    return fig


