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

    try:
        rep = int(name.split('_')[-2])
    except:
        raise ValueError(f"Can't get number of repeats from read '{name}'")
    return rep

def seq_generator(handle):
    """Generator that yields a sequence and its name from a FASTA file"""
    name, seq = '', ''
    for line in handle:
        if line.startswith('>'):
            if name != '':
                yield name, seq
            name = line.strip()
            seq = ''
        else:
            seq = seq + line.strip()
    if name != '':
        yield name, seq