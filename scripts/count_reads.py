#!/usr/env/bin python3

# count the nubmer of reads in particular file tyeps
# filetypes:
# fastq - one read per four lines
# fasta - count number of '^>' lines
# variant tsv - count distinct read ids
# pivoted variant tsv - count number of lines
# distinct read counts - count number of lines

import argparse
import os
from scripts.utils import use_open

FA_EXTS = {'.fa', '.fna', '.fasta'}
FQ_EXTS = {'.fq', '.fastq'}

def main():
    
    parser = argparse.ArgumentParser(description='Count the number of reads in files')
    parser.add_argument('--files', type=str, nargs='*', help='files to count reads in', default=[])
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()

    # count the number of reads in each file
    # generate lines for tsv with (filename, filetype, count)
    lines = []
    for file in args.files:
        print(f"Working on file {file} of type {get_file_type(file)}")
        file_type = get_file_type(file)
        if file_type == 'fasta':
            count = count_fasta(file)
        elif file_type == 'fastq':
            count = count_fastq(file)
        elif file_type == 'variant_tsv':
            count = count_variant_tsv(file)
        elif file_type == 'pivoted_tsv':
            count = count_pivoted_tsv(file)
        elif file_type == 'distinct_read_counts':
            count = count_pivoted_tsv(file)
        else:
            raise ValueError('Unknown file type: {}'.format(file_type))

        lines.append([file, file_type, count])

    if args.output is None:
        for line in lines:
            print('\t'.join(map(str, line)))
        return

    # write the lines to the output file
    with open(args.output, 'wt') as f:
        for line in lines:
            f.write('\t'.join(map(str, line)) + '\n')

def get_file_type(file):

    # check file extension
    ext = os.path.splitext(file)[1]
    if ext == '.gz':
        ext = os.path.splitext(os.path.splitext(file)[0])[1]
    if ext in FA_EXTS:
        return 'fasta'
    if ext in FQ_EXTS:
        return 'fastq'

    with use_open(file, 'rt') as f:
        try:
            for i, l in enumerate(f):
                parts = l.split('\t')
                # check the first line of the file to determine the file type
                if l.startswith('>'):
                    return 'fasta'
                if l.startswith('@'):
                    return 'fastq'
                # long form of variant file
                if l.startswith('reference_name'):
                    return 'variant_tsv'
                # short form of variant file
                if l.startswith('query_name'):
                    return 'variant_tsv'
                # pivoted reads file
                if l.startswith('read_id'):
                    return 'pivoted_tsv'
                # counts of distinct reads
                if len(parts) == 2:
                    return 'distinct_read_counts'
                if 'count' in parts:
                    return 'distinct_read_counts'
                # if we have read the first 100 lines and not found a file type, give up
                if i > 100:
                    break
        # if not a text file, return unknown
        except UnicodeDecodeError:
            pass

    # if we have not found a file type, return unknown  
    return 'unknown'


def count_fasta(file):

    count = 0
    # count lines that start with > - one for each sequence
    with use_open(file, 'rt') as f:
        try:
            for i, l in enumerate(f):
                if l.startswith('>'):
                    count += 1
        except UnicodeDecodeError:
            pass
    return count

def count_fastq(file):

    # count lines starting with @
    count = 0
    with use_open(file, 'rt') as f:
        try:
            for i, l in enumerate(f):
                if l != '' and l[0] == '@':
                    count += 1
        except UnicodeDecodeError:
            pass

    return count
    

def count_variant_tsv(file):
    
    prev_read = ''
    count = -1 # account for header
    # iterate over lines of file
    with use_open(file, 'rt') as f:
        for l in f:
            # get read name from line
            parts = l.split('\t')
            read = parts[2]
            # count every time we see a different read
            if read != prev_read:
                count += 1
                prev_read = read

    return count

def count_pivoted_tsv(file):

    count = -1 # account for header
    with use_open(file, 'rt') as f:
        try:
            for line in f:
                if line != '':
                    count += 1
        except UnicodeDecodeError:
            count = 0

    return count 

if __name__ == '__main__':
    main()
