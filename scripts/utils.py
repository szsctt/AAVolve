import gzip

def use_open(filename, *args, **kwargs):
    """
    Use gzip.open if file is gzipped
    """
    
    if filename.endswith('.gz'):
        return gzip.open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)
    
def get_repeats_from_r2c2_name(name):

    return int(name.split('_')[-2])