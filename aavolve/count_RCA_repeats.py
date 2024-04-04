#!/usr/bin/env python3

# count the number of repeats in a C3POa output file

import argparse
from sys import argv

from aavolve.utils import use_open, get_repeats_from_r2c2_name, seq_generator

def main(argv):
    
        args = get_args(argv)
    
        repeats = count_repeats(args.input)

        write_output(args.output, repeats)

def get_args(argv):
    parser = argparse.ArgumentParser(description='Count the number of repeats in a C3POa output file')
    parser.add_argument('--input', '-i', help='Input file', required=True)
    parser.add_argument('--output', '-o', help='Output file', required=True)
    return parser.parse_args(argv)

def count_repeats(infile):
     
    repeats = {}
    with use_open(infile, 'rt') as in_f:
        # iterate over sequences
        for name, _ in seq_generator(in_f):
            # get number of repeats from each read
            reps = get_repeats_from_r2c2_name(name)
            # increment count of repeats
            if reps in repeats:
                repeats[reps] += 1
            else:
                repeats[reps] = 1

    return repeats


def write_output(outfile, repeats):
    with use_open(outfile, 'wt') as out_f:
        # write header
        out_f.write('Repeats\tCount\n')
        for rep, count in repeats.items():
            out_f.write(f'{rep}\t{count}\n')

if __name__ == '__main__':
    main(argv[1:])