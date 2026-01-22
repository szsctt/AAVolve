import os
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import re
import hashlib

from aavolve.utils import MAX_SEQS

def display_warning_file(path, title="Warning"):
    if path is None or path == "":
        return
    if not os.path.exists(path):
        return
    with open(path, "rt") as handle:
        text = handle.read().strip()
    if not text:
        return
    try:
        from IPython.display import Markdown, display
    except Exception:
        print(f"{title}: {text}")
        return
    display(Markdown(f"> **{title}**  \n" + text.replace("\n", "  \n")))

#### counts of reads ####

def import_read_count_data(df_file, seq_type):

    df = pd.read_csv(df_file, delimiter='\t', header=None, names=['filename', 'file_type', 'Count'])

    df['File type'] = df.apply(lambda row: assign_file_type(row['filename'], row['file_type'], seq_type), axis=1)

    assert all(df['File type'].notnull())

    input_count = df[df['File type'] == 'Input']['Count'].tolist()
    assert len(input_count) == 1
    input_count = input_count[0]

    df['Fraction of reads'] = df['Count'] / input_count

    assert all(df['Fraction of reads'] <= 1.0)
   
    return df

def assign_file_type(file_name, file_type, seq_type):
    if file_type == "fastq":
        return 'Input'
    elif file_type == "fasta":
        if seq_type == "np-cc":
            dir = os.path.dirname(file_name)
            dir2 = os.path.dirname(dir)
            dir3 = os.path.dirname(dir2)
            if os.path.basename(dir3) == "c3poa":
                return "Consensus"
            elif os.path.basename(dir) == "c3poa_filt":
                return "Filtered by repeats"
        elif seq_type == "sg":
            return "Input"   
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
    raise ValueError(f"File type {file_type} not recognized")

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

    pass_frac = df[df['File type'] == "Filtered non-parental variants"]['Fraction of reads'].to_list()[0]

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

    if len(df) == 0:
        print("No reads passing all filters.")
        return None

    # reverse order of rows to get most frequent at top
    #df = df.iloc[::-1]

    # get counts
    counts = df['count']
    
    # get a list of unique parents
    # drop count column
    df = df.drop(columns=['count'])

    # ensure df is a DataFrame (drop may return Series in degenerate cases)
    if isinstance(df, pd.Series):
        df = df.to_frame().T

    # elementwise replace 'non_parental_\d+' and collapse comma-containing values to 'multiple'
    def norm_val(x):
        s = str(x)
        s = re.sub(r'non_parental_\d+', 'non parental', s)
        if ',' in s:
            return 'multiple'
        return s

    df = df.applymap(norm_val)

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
            showscale=False,
        )

    p2 = go.Scatter(x=counts, 
                    y=df_nums.index, 
                    mode='lines+markers', 
                    line_color='black', 
                    marker=dict(color='black'), 
                    showlegend=True,
                    name="Count",
                    legendgroup="Metrics",
                    legendgrouptitle_text="Metrics",
        )
    fig = make_subplots(rows=1, cols=6, specs = [[{"colspan": 1}, {"colspan": 5}, None, None, None, None]], shared_yaxes=True, horizontal_spacing=0.05, vertical_spacing=0.05)
    fig.add_trace(p2, row=1, col=1)
    fig.add_trace(p1, row=1, col=2)

    # Add a single combined legend for parent colors (instead of a separate colorbar).
    for parent in parents:
        fig.add_trace(
            go.Scatter(
                x=[0],
                y=[0],
                mode="markers",
                marker=dict(
                    size=10,
                    color=color_dict[parent],
                    line=dict(width=0.5, color="black"),
                ),
                name=str(parent),
                visible="legendonly",
                legendgroup="Parents",
                legendgrouptitle_text="Parents",
            ),
            row=1,
            col=2,
        )
    fig['layout']['xaxis']['title'] = 'Count'
    fig['layout']['yaxis']['title'] = 'Read'
    fig['layout']['xaxis2']['title'] = 'Position in reference'

    # Improve readability and prevent legend clipping in reports.
    max_parent_len = max((len(p) for p in parents), default=0)
    right_margin = min(520, 180 + int(max_parent_len * 6.5))
    fig.update_layout(
        margin=dict(l=60, r=right_margin, t=40, b=90),
        legend=dict(x=1.02, xanchor="left", y=1, yanchor="top", font=dict(size=10)),
    )
    fig.update_xaxes(automargin=True)
    fig.update_xaxes(tickangle=90, automargin=True, nticks=30, row=1, col=2)

    return fig

def plot_breakpoints(breakpoints_file, counts_file, seq_type):

    # get total reads passing all filters
    counts = import_read_count_data(counts_file, seq_type)
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

    # get location of each variant in the reference
    df['variant'] = numeric_position(df['variant'])
    # sort by position
    df['variant_pos'] = df['variant'].str.extract(r'(^\d+)').astype(int)
    df = df.sort_values('variant_pos')

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
                  xaxis={'categoryorder':'array', 'categoryarray':df['variant']},
                 )

    # Improve readability and avoid legend clipping.
    fig.update_xaxes(tickangle=90, automargin=True, nticks=40)
    max_parent_len = max((len(p) for p in parent_colors_dict.keys()), default=0)
    right_margin = min(420, 120 + int(max_parent_len * 6.5))
    fig.update_layout(
        margin=dict(l=60, r=right_margin, t=60, b=140),
        legend=dict(x=1.02, xanchor='left', y=1, yanchor='top', font=dict(size=10)),
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

    # Get unique parents in file order (not sorted).
    parents = []
    seen = set()
    for parent in df['parent'].tolist():
        if parent in seen:
            continue
        seen.add(parent)
        parents.append(parent)

    # Ensure the special categories are present and placed at the end.
    for special in ("non parental", "multiple"):
        if special in parents:
            parents.remove(special)
    parents.extend(["non parental", "multiple"])

    # Choose a palette based on how many distinct parents we need to show.
    #
    # Keep this stable across Plotly versions by using an explicit Prism palette
    # (some Plotly versions expose Prism as RGB strings and with a shorter length).
    prism_hex = [
        '#FD3216', '#00FE35', '#6A76FC', '#FED4C4', '#FE00CE', '#0DF9FF', '#F6F926', '#FF9616',
        '#479B55', '#EEA6FB', '#DC587D', '#D626FF', '#6E899C', '#00B5F7', '#B68E00', '#C9FBE5',
    ]
    total = len(parents)
    if total <= len(px.colors.qualitative.Plotly):
        colors = px.colors.qualitative.Plotly[:total]
    elif total <= len(prism_hex):
        colors = prism_hex[:total]
    else:
        # Fall back to a continuous colorscale sampled across [0, 1].
        xs = [i / (total - 1) for i in range(total)]
        colors = px.colors.sample_colorscale(px.colors.sequential.Turbo, xs)

    return dict(zip(parents, colors))

def numeric_position(col):

    col = col.astype(str).str.replace(":sub", "", regex=True)
    col = col.str.replace(":ins", " ", regex=True)
    col = col.str.replace(":del", "  ", regex=True)

    return col
