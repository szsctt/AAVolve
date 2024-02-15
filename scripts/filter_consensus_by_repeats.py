#!/usr/bin/env python3

# Filter C3POa output by number of repeats

import argparse
from sys import argv

from scripts.utils import use_open, get_repeats_from_r2c2_name

def main(argv):

    parser = argparse.ArgumentParser(description='Filter C3POa output by number of repeats')
    parser.add_argument('--input', '-i', help='Input file', required=True)
    parser.add_argument('--output', '-o', help='Output file', required=True)
    parser.add_argument('--min-repeats', '-m', type=int, default=2, help='Minimum number of repeats')
    args = parser.parse_args(argv)

    filter(args.input, args.output, args.min_repeats)


def filter(infile, outfile, min_repeats):
    """Filter C3POa output by number of repeats"""

    if min_repeats < 1:
        raise ValueError('min_repeats must be >= 1')

    with use_open(infile, 'rt') as in_f, use_open(outfile, 'wt') as out_f:
        for name, seq in seq_generator(in_f):
            # check number of repeats
            repeats = int(name.split('_')[-2])
            # write to output file if repeats >= min_repeats
            if repeats >= min_repeats:
               out_f.write(f'{name}\n{seq}\n')

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
        

if __name__ == '__main__':
    main(argv[1:])

