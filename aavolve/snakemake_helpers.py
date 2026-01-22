import numpy as np
from snakemake.io import expand
import hashlib
from pathlib import Path
import re

# general helpers
def get_column_by_sample(wildcards, samples, column_name):

    assert len(samples.sample_name) == len(samples.sample_name.unique()), "Sample names are not unique"
    return {k:v for k, v in zip(samples.sample_name, samples[column_name])}[wildcards.sample]

def get_column_by_parent(wildcards, samples, column_name):
    # check that column exists
    if column_name not in samples.columns:
        raise KeyError(f"Column '{column_name}' not found in samples DataFrame")
    # check that each parent_name corresponds to the same column value
    test = samples[['parent_name', column_name]].drop_duplicates() # unique combinations of parent_name and column_name
    counts = test.groupby('parent_name').size().reset_index(name='count') # each parent should only appear once
    assert all(counts['count'] == 1)
    return {k:v for k, v in zip(samples.parent_name, samples[column_name])}[wildcards.sample]


def is_fastq(file):
    return any((file.endswith('.fastq'), file.endswith('.fastq.gz'), file.endswith('.fq'), file.endswith('.fq.gz')))


def get_anchor_length(wildcards, samples):
    """Return anchor length configured for the sample/parent."""

    column = 'anchors'
    if column not in samples.columns:
        return 0

    if wildcards.sample in set(samples.sample_name):
        value = get_column_by_sample(wildcards, samples, column)
    elif wildcards.sample in set(samples.parent_name):
        value = get_column_by_parent(wildcards, samples, column)
    else:
        raise ValueError(f"Sample {wildcards.sample} not found in samples DataFrame")

    if value is None:
        return 0

    if isinstance(value, (float, np.floating)) and np.isnan(value):  # type: ignore[arg-type]
        return 0

    if isinstance(value, bool):
        raise ValueError(
            f"Anchor length must be a non-negative integer. Found boolean {value} for sample {wildcards.sample}"
        )

    try:
        length = int(value)
    except (TypeError, ValueError) as err:
        raise ValueError(
            f"Anchor length must be a non-negative integer. Found value {value!r} for sample {wildcards.sample}"
        ) from err

    if length < 0:
        raise ValueError(f"Anchor length must be a non-negative integer. Found {length} for sample {wildcards.sample}")

    return length


def get_anchor_seed(wildcards, samples):
    """
    Generate a consistent seed for anchor generation based on parent+reference+anchor_length.
    Samples with the same parent, reference, and anchor length will share the same anchors.
    """
    # Get parent name
    if wildcards.sample in set(samples.parent_name):
        parent_name = wildcards.sample
    else:
        parent_name = get_column_by_sample(wildcards, samples, 'parent_name')
    
    # Get reference name
    if wildcards.sample in set(samples.parent_name):
        reference_name = get_column_by_parent(wildcards, samples, 'reference_name')
    else:
        reference_name = get_column_by_sample(wildcards, samples, 'reference_name')
    
    # Get anchor length
    anchor_length = get_anchor_length(wildcards, samples)
    
    # Create a consistent seed from parent+reference+length
    seed = f"{parent_name}_{reference_name}_{anchor_length}"
    return seed


def anchors_enabled(wildcards, samples):
    return get_anchor_length(wildcards, samples) > 0


def input_pair_id(parent_file: str, reference_file: str) -> str:
    """Return a stable identifier for a (parent_file, reference_file) pair.

    This is used to validate inputs once per unique file combination.
    """

    def _basename_without_compression(path: str) -> str:
        name = Path(path).name
        return name[:-3] if name.endswith(".gz") else name

    def _stem(path: str) -> str:
        return Path(_basename_without_compression(path)).stem

    def _slug(text: str) -> str:
        slug = re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-")
        return slug or "group"

    parent_label = _slug(_stem(parent_file))
    reference_label = _slug(_stem(reference_file))

    # Add a short, stable hash suffix to avoid collisions (e.g. same basename in different dirs).
    key = f"{parent_file}|{reference_file}"
    suffix = hashlib.md5(key.encode("utf-8")).hexdigest()[:8]
    return f"{parent_label}__{reference_label}__{suffix}"


def build_input_validation_targets(samples_df):
    """Build unique input-pair maps and the corresponding validation targets.

    Returns:
        input_validation_map: dict[pair_id] -> (parent_file, reference_file)
        input_validation_samples: dict[pair_id] -> list[sample_name]
        input_validation_targets: list of dummy files to add to rule all

    Notes:
        This function does not read any of the referenced files.
    """

    input_validation_map = {}
    input_validation_samples = {}

    for parent_file, reference_file, sample_name in zip(
        samples_df.parent_file, samples_df.reference_file, samples_df.sample_name
    ):
        pair_id = input_pair_id(parent_file, reference_file)
        input_validation_map[pair_id] = (parent_file, reference_file)
        input_validation_samples.setdefault(pair_id, []).append(sample_name)

    input_validation_targets = [
        f"out/qc/input-checks/{pair_id}.txt" for pair_id in sorted(input_validation_map)
    ]

    return input_validation_map, input_validation_samples, input_validation_targets


def build_group_report_targets(samples_df):
    """Return group report targets keyed by unique parent/reference pairs.

    Notes:
        This function does not read any of the referenced files.
    """

    input_validation_map, input_validation_samples, _ = build_input_validation_targets(samples_df)
    group_report_targets = [
        f"out/reports/group_reports/{pair_id}_report.html" for pair_id in sorted(input_validation_map)
    ]
    return input_validation_map, input_validation_samples, group_report_targets


def build_non_parental_variant_group_maps(samples_df, input_validation_samples):
    """Return maps needed to share non-parental variants within parent/reference groups.

    This reuses the group identifiers already computed for input validation / group reports
    by consuming the ``input_validation_samples`` mapping produced by
    ``build_input_validation_targets`` instead of recomputing pair IDs.

    Args:
        samples_df: samples DataFrame (must include columns sample_name, include_non_parental).
        input_validation_samples: dict[pair_id] -> list[sample_name] for all samples in each group
            (output of ``build_input_validation_targets``).

    Returns:
        sample_to_pair_id: dict mapping sample_name -> pair_id for all samples
        pair_id_to_samples: dict mapping pair_id -> list[sample_name] for samples with include_non_parental=True
    """

    include_non_parental_by_sample = {
        str(sample_name): bool(include_non_parental)
        for sample_name, include_non_parental in zip(samples_df.sample_name, samples_df.include_non_parental)
    }

    sample_to_pair_id = {}
    pair_id_to_samples = {}

    for pair_id, sample_names in input_validation_samples.items():
        for sample_name in sample_names:
            sample_name = str(sample_name)
            sample_to_pair_id[sample_name] = pair_id
            if include_non_parental_by_sample.get(sample_name, False):
                pair_id_to_samples.setdefault(pair_id, []).append(sample_name)

    return sample_to_pair_id, pair_id_to_samples


def get_non_parental_high_freq_variants_for_sample(
    wildcards,
    samples_df,
    sample_to_pair_id,
    non_parental_variant_groups,
):
    """Return the appropriate non-parental variant file for a sample.

    If a sample has include_non_parental=True and shares (parent_file, reference_file)
    with other include_non_parental samples, return the combined non-parental variants
    across the group; otherwise return the per-sample file.
    """

    include_non_parental = bool(
        samples_df.loc[samples_df.sample_name == wildcards.sample, "include_non_parental"].iloc[0]
    )
    if not include_non_parental:
        return f"out/variants/frequency/{wildcards.sample}_high.tsv.gz"

    pair_id = sample_to_pair_id[wildcards.sample]
    group_samples = non_parental_variant_groups.get(pair_id, [])
    if len(group_samples) <= 1:
        return f"out/variants/frequency/{wildcards.sample}_high.tsv.gz"

    return f"out/variants/frequency/groups/{pair_id}_high.tsv.gz"


def get_anchor_sequences_path(wildcards):
    return f"out/anchors/{wildcards.sample}.fasta"


def get_reads_for_anchor_input(wildcards, samples):
    if wildcards.sample in set(samples.parent_name):
        return get_reads(wildcards, samples)
    if trim_is_enabled(wildcards, samples):
        return f"out/trimmed/{wildcards.sample}.trimmed.gz"
    return get_reads(wildcards, samples)


def get_anchor_reads_suffix(wildcards, samples):
    base_reads = get_reads(wildcards, samples)

    if wildcards.sample in set(samples.parent_name):
        return '.fasta.gz' if base_reads.endswith('.gz') else '.fasta'

    if trim_is_enabled(wildcards, samples):
        return '.fastq.gz'

    if is_fastq(base_reads):
        return '.fastq.gz' if base_reads.endswith('.gz') else '.fastq'

    return '.fasta.gz' if base_reads.endswith('.gz') else '.fasta'


def get_anchored_reads_path(wildcards, samples):
    suffix = get_anchor_reads_suffix(wildcards, samples)
    return f"out/anchors/reads/{wildcards.sample}{suffix}"


def get_anchored_reference_path(wildcards):
    return f"out/anchors/references/{wildcards.sample}.fasta"


def _get_reference_base(wildcards, samples):
    if wildcards.sample not in set(samples.sample_name) | set(samples.parent_name):
        raise ValueError(f"Sample {wildcards.sample} not found")

    parents = {}
    for k, v in zip(samples['parent_name'], samples['reference_file']):
        parents[k] = v

    if wildcards.sample in parents.keys():
        return parents[wildcards.sample]

    return get_column_by_sample(wildcards, samples, 'reference_file')


def minimap2_params_with_default(wildcards, samples):
    default = "-x map-hifi -B 1.5 --end-bonus 5"
    if wildcards.sample not in set(samples.sample_name) | set(samples.parent_name):
        raise ValueError(f"Sample {wildcards.sample} not found in samples DataFrame")
    if wildcards.sample in set(samples.sample_name) and wildcards.sample in set(samples.parent_name):
        raise ValueError(f"Sample {wildcards.sample} is both a sample and a parent, please clarify which one to use")
    
    # get params from samples
    user_params = default
    if wildcards.sample in set(samples.parent_name):
        # If the sample is a parent, use the parent's minimap2 parameters
        try:
            user_params = get_column_by_parent(wildcards, samples, "minimap2_params")
        except KeyError:
            pass
    if wildcards.sample in set(samples.sample_name):
        # If the sample is not a parent, use the sample's minimap2 parameters
        try:
            user_params = get_column_by_sample(wildcards, samples, "minimap2_params")
        except KeyError:
            pass
    
    # add default parameters if not provided
    if not isinstance(user_params, str):
        raise ValueError(f"Expected 'minimap2_params' to be a string, got {type(user_params)}")
    # Ensure -x or --preset is present
    if "-x" not in user_params and "--preset" not in user_params:
        user_params += f" -x map-hifi"
    # Ensure -B is present
    if "-B" not in user_params and "--score-N" not in user_params:
        user_params += " -B 1.5"
    # Ensure --end-bonus is present
    if "--end-bonus" not in user_params:
        user_params += " --end-bonus 5"
    return user_params.strip()

#### align ####

def get_reads(wildcards, samples):
    """
    Get appropriate reads for wildcards.sample
    Either parental references,
    output from C3POa for nanopore R2C2 reads
    or just fastq otherwise
    """
    # make a dictionary of parents
    parents = {}
    for k, v in zip(samples['parent_name'], samples['parent_file']):
        parents[k] = v
    
    # if one of the parents, return parent sequences
    if wildcards.sample in parents.keys():
        return parents[wildcards.sample]
    
    # get sequencing technology
    tech = get_column_by_sample(wildcards, samples, 'seq_tech')
   
    # if nanopore r2c2, return consensus reads
    if tech == 'np-cc':
        
        # check if we want to filter for repeats
        repeats = get_column_by_sample(wildcards, samples, 'min_reps')

        if np.isnan(repeats):
            return f"out/c3poa/{wildcards.sample}/splint/R2C2_Consensus.fasta.gz"

        else:
            return f"out/c3poa_filt/{wildcards.sample}.fasta.gz"

    # otherwise, just return reads
    return get_column_by_sample(wildcards, samples, 'read_file')

def get_reference(wildcards, samples):
    """Return the raw reference path (without anchors)."""

    return _get_reference_base(wildcards, samples)


def get_reference_for_align(wildcards, samples):
    """Return reference path to use during alignment and variant extraction."""

    if anchors_enabled(wildcards, samples):
        return get_anchored_reference_path(wildcards)

    return _get_reference_base(wildcards, samples)


def trim_is_enabled(wildcards, samples):
    """Return True if adapter trimming is requested for the sample."""

    # parents never undergo trimming
    if wildcards.sample in set(samples.parent_name):
        return False

    try:
        trim_value = get_column_by_sample(wildcards, samples, 'trim')
    except KeyError:
        return False

    if not isinstance(trim_value, bool):
        raise ValueError(
            f"Unexpected value for 'trim' column for sample '{wildcards.sample}': {trim_value!r}. "
            "Expected boolean True/False."
        )

    return bool(trim_value)


def get_reads_for_align(wildcards, samples):
    """Return the appropriate input for minimap2 alignment, considering trimming."""

    if anchors_enabled(wildcards, samples):
        return get_anchored_reads_path(wildcards, samples)

    reads = get_reads(wildcards, samples)

    if not trim_is_enabled(wildcards, samples):
        return reads

    return f"out/trimmed/{wildcards.sample}.trimmed.gz"


def get_linked_adapters_for_sample(wildcards, samples):
    """Construct the linked-adapter specification for cutadapt, validating inputs."""

    if not trim_is_enabled(wildcards, samples):
        return None

    try:
        adapter_5 = get_column_by_sample(wildcards, samples, 'adapter_5')
        adapter_3 = get_column_by_sample(wildcards, samples, 'adapter_3')
    except KeyError as err:
        missing = err.args[0]
        raise ValueError(
            f"Trimming is enabled for sample '{wildcards.sample}' but column '{missing}' is missing. "
            "Please ensure both 'adapter_5' and 'adapter_3' are present in the samples configuration."
        ) from err

    if not isinstance(adapter_5, str):
        raise ValueError(
            f"Trimming is enabled for sample '{wildcards.sample}', but adapter_5 is not a string (value: {adapter_5!r})."
        )
    if not isinstance(adapter_3, str):
        raise ValueError(
            f"Trimming is enabled for sample '{wildcards.sample}', but adapter_3 is not a string (value: {adapter_3!r})."
        )

    adapter_5 = adapter_5.strip()
    adapter_3 = adapter_3.strip()

    if adapter_5 == '' or adapter_3 == '':
        raise ValueError(
            f"Trimming is enabled for sample '{wildcards.sample}', but adapter sequences are missing or empty. "
            "Populate 'adapter_5' and 'adapter_3' in the samples configuration."
        )

    return f"{adapter_5}...{adapter_3}"

#### variants ####

# rules num_parents
def get_parents(wildcards, samples):
    
    if wildcards.sample not in set(samples.parent_name):
        raise ValueError(f"Parent {wildcards.sample} not found")

    parents = {k:v for k,v in zip(samples.parent_name, samples.parent_file)}[wildcards.sample]
    
    # check there is a value for the parent file
    if parents is None:
        raise ValueError(f"Sample {wildcards.sample} does not have a parent")
    try:
        if np.isnan(parents):
            raise ValueError(f"Sample {wildcards.sample} does not have a parent")
    except TypeError:
        pass

    return parents


def fill_parents(wildcards, samples, filename):
    parent_name = get_column_by_sample(wildcards, samples, "parent_name")
    return expand(filename, sample=parent_name)

#### tranform_variants ###

# rule dmat
def get_dmat_input(wildcards, nt_seq_file, aa_seq_file):
    if wildcards.seq_type == "nt-seq":
        return nt_seq_file
    elif wildcards.seq_type == "aa-seq":
        return aa_seq_file
    else:
        raise ValueError("seq_type must be nt-seq or aa-seq")

# rule count_reads
def get_reads_for_counting(wildcards, samples, consensus, con_filt):
    """
    Get the appropriate files for counting original reads
    If not R2C2 data, just return result of get_reads(wildcards)

    If R2C2 data, return original reads and consensus reads
    If we filtered, also return filtered reads
    """

    # get input reads
    reads = [{k:v for k, v in zip(samples.sample_name, samples.read_file)}[wildcards.sample]]
    
    # get sequencing technology
    tech = {k:v for k, v in zip(samples.sample_name, samples.seq_tech)}[wildcards.sample]
    # if not R2C2, just return reads
    if tech != 'np-cc':
        return reads

    # for R2C2 add consensus reads
    reads.append(consensus)
    # and filtered consensus reads if they exist
    filt = {k:v for k, v in zip(samples.sample_name, samples.min_reps)}[wildcards.sample]
    if filt is not None and np.isnan(filt) == False:
        reads.append(con_filt)
    
    return reads

# rule count_reads
def format_input_reads(input):
    """
    Format input reads for counting
    """
    fastqs = []
    fastas = []
    
    # check if input was fastq (most reads) or fasta (sanger sequencing)
    for i in input:
        if is_fastq(i):
            fastqs.append(i)
        else:
            fastas.append(i)
    
    # format output
    args = ''
    if len(fastas) > 0:
        args = args + f"--fasta-files {' '.join(fastas)} "
    if len(fastqs) > 0:
        args = args + f"--fastq-files {' '.join(fastqs)}"

    return args
