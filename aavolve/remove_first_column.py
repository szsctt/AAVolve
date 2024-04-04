#!/usr/bin/env python3

# remove the first column from a file

import argparse

from aavolve.utils import use_open

def main():
    
    parser = argparse.ArgumentParser(description='Remove the first column from a file')
    parser.add_argument('-i', '--input', help='input file', required=True)
    args = parser.parse_args()

    with use_open(args.input, 'rt') as f:

        # skip header
        next(f)

        for line in f:
                print(remove_first_column(line), end='')

def remove_first_column(line):
    
    cols = line.split('\t')

    if len(cols) < 2:
        return ''

    return '\t'.join(cols[1:])


if __name__ == "__main__":
    main()