#!/usr/bin/env python3

# get the first and last variant from a file with variants from parents

import argparse
from sys import argv

from scripts.utils import read_variant_file

def main(args):

    arg = get_args(args)

    first, last = get_first_last_variant(arg.input)

    write_output(arg.output, first, last)


def get_first_last_variant(filename):

    first = None
    last = None

    for variant in read_variant_file(filename):
        pos = int(variant['pos']) if '_' not in variant['pos'] else int(variant['pos'].split('_')[0])
        if first is None or pos < first:
            first = pos
        if last is None or pos > last:
            last = pos

    assert first is not None
    assert last is not None

    return first, last


def write_output(filename, first, last):

    with open(filename, 'w') as f:
        f.write(str(first) + '\n')
        f.write(str(last) + '\n')

def get_args(args):

    parser = argparse.ArgumentParser(description="Get the first and last variant from a file with variants from parents")
    parser.add_argument('--input', '-i', help="File with variants from parents", required=True)
    parser.add_argument('--output', '-o', help="Output file with first and last variant", required=True)
    return parser.parse_args(args)

if __name__ == "__main__":

    main(argv[1:])