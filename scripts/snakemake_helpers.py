def get_column_by_sample(wildcards, samples, column_name):

    assert len(samples.sample_name) == len(samples.sample_name.unique()), "Sample names are not unique"
    return {k:v for k, v in zip(samples.sample_name, samples[column_name])}[wildcards.sample]