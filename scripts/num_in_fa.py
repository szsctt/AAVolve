#!/usr/bin/env python3

# count the number of records in a fasta file

import argparse
from sys import argv

from scripts.utils import seq_generator

def main(args):
    
        arg = get_args(args)
    
        n = count_records(arg.input)

        write_output(arg.output, n)
    

def count_records(in_file):
      
    count = 0
    with open(in_file, 'r') as f:
        for _ in seq_generator(f):
            count += 1

    return count

def write_output(out_file, n):
    with open(out_file, 'w') as f:
        f.write(str(n) + '\n')


def get_args(args):
        
    parser = argparse.ArgumentParser(description="Count the number of records in a fasta file")
    parser.add_argument('--input', '-i', help="Input fasta file", required=True)
    parser.add_argument('--output', '-o', help="Output file with the number of records", required=True)
    return parser.parse_args(args)
    

if __name__ == "__main__":
    main(argv[1:])