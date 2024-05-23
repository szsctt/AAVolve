import os
import numpy as np
import pandas as pd

PARENTDIR = 'out/references/parents'
os.makedirs(PARENTDIR, exist_ok=True)
DEFAULT_MINREPS = 3
REQUIRED_COLUMNS =  ('sample_name', 'parent_name', 'reference_name', 'seq_tech', 'read_file', 'parent_file', 'reference_file')
SEQ_TECHS = ['np', 'np-cc', 'pb', 'pb-hifi', 'sg']
DEFAULT_FREQ = 0.2
DEFAULT_INCLUDE_NON = False
DEFAULT_GROUP_VARS = False

def get_name(filename):
    return os.path.splitext(os.path.basename(filename))[0]

def get_first_parent(filename):
    """
    Get first parent from fasta file and save in separate file
    """
    
    with open(filename, 'r') as f:
        # read first line
        first_line = f.readline()
        if first_line[0] != '>':
            raise Exception(f"First line of {filename} does not start with '>'")
        else:
            parent_name = first_line[1:].strip().split()[0]
            parent_file = os.path.join(PARENTDIR, parent_name + '.fa')
            with open(parent_file, 'w') as p:
                p.write(first_line)
                for line in f:
                    if line[0] == '>':
                        break
                    p.write(line)

    return parent_name, parent_file

def get_command_options(config):
    """
    Get options from the config file
    """

    samples = dict()
    # required options
    required = {'read_file', 'parent_file'}
    for r in required:
        if r not in config:
            raise Exception(f"Please include the command line argument --config {r}=<path to {r}>")
        else:
            samples[r] = config[r]

    # deal with reference - if not specified, use the first parent from the parent file
    if 'reference_file' not in config:
        if not os.path.isfile(config['parent_file']):
            raise Exception(f"Parent file does not exist: {config['parent_file']}")
        
        # get first parent from parent file
        ref_name, ref_file = get_first_parent(config['parent_file'])
        print(f"Reference not specified: using the first parent ('{ref_name}') from {config['parent_file']}")
        samples['reference_file'] = ref_file
        config['reference_file'] = ref_file
        
        # if name not specified, use name of first parent 
        if 'refrence_name' not in config:
            samples['reference_name'] = ref_name
            config['reference_name'] = ref_name
        else:
            samples['reference_name'] = config['reference_name']
    # if reference is specified
    else:
        samples['reference_file'] = config['reference_file']
    
        if 'reference_name' not in config:
            samples['reference_name'] = get_name(config['reference_file'])
            config['reference_name'] = get_name(config['reference_file'])
        else:
            samples['reference_name'] = config['reference_name']
    
    # sequencing technology
    if 'seq_tech' not in config:
        raise Exception("Please include the command line argument --config seq_tech=<np, np-cc or pb>")
    else:
        samples['seq_tech'] = config['seq_tech']

    # for seq_tech = np-cc, need to specify splint file
    if 'splint_file' not in config:
        config['splint_file'] = None
    if config['seq_tech'] == 'np-cc' and config['splint_file'] is None:
        raise Exception("Please include the command line argument --config splint_file=<path to splint file>")
    else:
        samples['splint_file'] = config['splint_file']

    # optional options
    optional = {'sample_name':'read_file', 'parent_name':'parent_file', 'reference_name':'reference_file'}
    for col, rep in optional.items():
        if col not in config:
            samples[col] = get_name(config[rep])
        else:
            samples[col] = config[col]
    
    if 'min_reps' in config:
        samples['min_reps'] = config['min_reps']


    return pd.DataFrame(samples, index=[0])

#### sanity checks on data ####
def check_data(samples):

    # required columns: sample_name, parent_name, reference_name, seq_tech, read_file, parent_file, reference_file
    for col in REQUIRED_COLUMNS:
        # check that all required columns are present
        if col not in samples.columns:
            raise Exception(f"Missing required column: {col}")
        # check that all required columns are non-empty
        if samples[col].isnull().any():
            raise Exception(f"Missing value in column: {col}")

    # check that all sample files exist
    for read_file in set(samples['read_file']):
        if not os.path.isfile(read_file):
            raise Exception(f"Sample file does not exist: {read_file}")

    # check that all samples have a unique name
    if len(samples['sample_name'].unique()) != len(samples):
        raise Exception("Sample names (column 'sample_name') must be unique")

    # check that all parent files exist
    for parent_file in set(samples['parent_file']):
        if not os.path.isfile(parent_file):
            raise Exception(f"Parent file does not exist: {parent_file}")

    # check that the same parent name always correpsonds to the same parent file
    if len(samples.groupby(['parent_name', 'parent_file'])) != len(samples.groupby('parent_name')):
        raise Exception("Each parent name (column 'parent_name') must always correspond to the same parent file (column 'parent_file')")
    if len(samples.groupby(['parent_name', 'parent_file'])) != len(samples.groupby('parent_file')):
        raise Exception("Each parent name (column 'parent_name') must always correspond to the same parent file (column 'parent_file')")
    
    # check that none of the parent names are the same as the sample names
    if len(set(samples['sample_name']).intersection(set(samples['parent_name']))) > 0:
        raise Exception("Sample names (column 'sample_name') must be different from parent names (column 'parent_name')")

    # check that all reference files exist
    for reference_file in set(samples['reference_file']):
        if not os.path.isfile(reference_file):
            raise Exception(f"Reference file does not exist: {reference_file}")

    #check that the same reference name always corresponds to the same parent file
    if len(samples.groupby(['reference_name', 'reference_file'])) != len(samples.groupby('reference_name')):
        raise Exception("Each reference name (column 'reference_name') must always correspond to the same reference file (column 'reference_file')")
    if len(samples.groupby(['reference_name', 'reference_file'])) != len(samples.groupby('reference_file')):
        raise Exception("Each reference name (column 'reference_name') must always correspond to the same reference file (column 'reference_file')")
    
    # check that none of the reference names are the same as the sample names
    if len(set(samples['sample_name']).intersection(set(samples['reference_name']))) > 0:
        raise Exception("Sample names (column 'sample_name') must be different from reference names (column 'reference_name')")

    # check that sequencing technology is either 'np', 'np-cc', 'pb', 'pb-hifi' or 'sg'
    if not set(samples['seq_tech']).issubset(set(SEQ_TECHS)):
        techs = [f"'{s}'" for s in SEQ_TECHS]
        techs = ', '.join(techs)
        raise Exception(f"Sequencing technology (column 'seq_tech') must be one of {techs}")
    
    #for 'nanopore-cc', must specify splint file
    for i, row in samples.iterrows():
        if row['seq_tech'] == 'np-cc':
            if 'splint_file' not in row or pd.isnull(row['splint_file']):
                raise Exception("For sequencing technology 'np-cc', must specify splint file (column 'splint_file')")
            elif not os.path.isfile(row['splint_file']):
                raise Exception(f"Splint file does not exist: {row['splint_file']}")


    # for nanopore-cc, if minimum reps isn't specified, set to default value
    if 'min_reps' not in samples.columns:
        samples['min_reps'] = [None]*len(samples)
    for i, row in samples.iterrows():
        if row['seq_tech'] == 'np-cc':

            # check that value is integer or null
            try:
                pd.to_numeric(row['min_reps'], errors='raise')
            except ValueError:
                raise ValueError(f"Minimum reps (column 'min_reps') must be an integer and at least 0: found value {row['min_reps']} in row {i}")

            # if null, fill with default value
            if row['min_reps'] is None or pd.isnull(row['min_reps']):
                samples.loc[i, 'min_reps'] = DEFAULT_MINREPS

            # check that value is integer
            elif row['min_reps'] != int(row['min_reps']):
                raise ValueError(f"Minimum reps (column 'min_reps') must be an integer and at least 0: found value {row['min_reps']} in row {i}")

            # check not less than 0
            elif row['min_reps'] < 0:
                raise ValueError(f"Minimum reps (column 'min_reps') must be an integer and at least 0: found value {row['min_reps']} in row {i}")
        else:
            samples.loc[i, 'min_reps'] = None


    # check if non_parental_freq is specified - otherwise fill with default
    if 'non_parental_freq' not in samples.columns:
        samples['non_parental_freq'] = [DEFAULT_FREQ]*len(samples) 

    # check that all non_parental_freq is set and between 0 and 1
    for i, row in samples.iterrows():
        
        # check that value is numeric or null
        try:
            pd.to_numeric(row['non_parental_freq'], errors='raise')
        except ValueError:
            raise ValueError(f"Column 'non_parental_freq' must be numeric and between 0 and 1: found value {row['non_parental_freq']} in row {i}")

        # not null, check value is between 0 and 1
        if (row['non_parental_freq'] is not None or not pd.isnull(row['non_parental_freq'])) and (row['non_parental_freq'] < 0 or row['non_parental_freq'] > 1):
            raise Exception(f"Column 'non_parental_freq' must be numeric and between 0 and 1: found value {row['non_parental_freq']} in row {i}")
        
        # if null, set to default
        elif row['non_parental_freq'] is None or pd.isnull(row['non_parental_freq']):
            samples.loc[i, 'non_parental_freq'] = DEFAULT_FREQ


    # check if include_freqs is specified - otherwise fill with default
    if 'include_non_parental' not in samples.columns:
        samples['include_non_parental'] = [DEFAULT_INCLUDE_NON]*len(samples)
    
    # check that all include_non_parental values are True or False
    for i, row in samples.iterrows():

        # check that values are boolean or null
        if not (row['include_non_parental'] is True or row['include_non_parental'] is False or
                row['include_non_parental'] is None or pd.isnull(row['include_non_parental'])):
                raise ValueError(f"Column 'include_non_parental' must be True, False or omitted: found value {row['include_non_parental']} in row {i}")
        
        # if null, fill with default
        if row['include_non_parental'] is None or pd.isnull(row['include_non_parental']):
            samples.loc[i, 'include_non_parental'] = DEFAULT_INCLUDE_NON

    # check if group_vars is specified - otherwise fill with default
    if 'group_vars' not in samples.columns:
        samples['group_vars'] = [DEFAULT_GROUP_VARS]*len(samples)

    # check that all group_vars values are True or False
    for i, row in samples.iterrows():

        # check that values are boolean or null
        if not (row['group_vars'] is True or row['group_vars'] is False or
                row['group_vars'] is None or pd.isnull(row['group_vars'])):
                raise ValueError(f"Column 'group_vars' must be True, False or omitted: found value {row['group_vars']} in row {i}")
        
        # if null, fill with default
        if row['group_vars'] is None or pd.isnull(row['group_vars']):
            samples.loc[i, 'group_vars'] = DEFAULT_GROUP_VARS

    # can't do both group_vars and include_non_parental for the same sample
    for i, row in samples.iterrows():

        if row['group_vars'] is True and row['include_non_parental'] is True:
            raise ValueError(f"Error in row {i}: Can't set both 'group_vars' and 'include_non_parental' to True")

    return samples

def get_samples(config):

    # read sample csv file
    if 'samples' not in config:
        samples = get_command_options(config)
    else:   
        samples = pd.read_csv(config['samples'])
        
    # do checks
    samples = check_data(samples)

    return(samples)
