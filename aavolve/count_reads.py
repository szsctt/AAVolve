#!/usr/env/bin python3

# count the nubmer of reads in particular file tyeps
# filetypes:
# fastq - one read per four lines
# fasta - count number of '^>' lines
# variant tsv - count distinct read ids
# pivoted variant tsv - count number of lines
# distinct read counts - count number of lines

import argparse
from aavolve.utils import use_open

FA_EXTS = {'.fa', '.fna', '.fasta'}
FQ_EXTS = {'.fq', '.fastq'}

def main():
    
    parser = argparse.ArgumentParser(description='Count the number of reads in files')
    parser.add_argument('--fastq-files', type=str, nargs='*', help='fastq files to count reads in', default=[])
    parser.add_argument('--fasta-files', type=str, nargs='*', help='fasta files to count reads in', default=[])
    parser.add_argument('--variant-read-ids', type=str, nargs='*', help='variant read id files to count reads in', default=[])
    parser.add_argument('--pivoted-tsv-files', type=str, nargs='*', help='pivoted tsv files to count reads in', default=[])
    parser.add_argument('--distinct-read-counts-files', type=str, nargs='*', help='distinct read counts files to count reads in', default=[])
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()

    # count the number of reads in each file
    # generate lines for tsv with (filename, filetype, count)
    lines = []
    for file in args.fastq_files:
        count = count_fastq(file)
        lines.append([file, 'fastq', count])
    for file in args.fasta_files:
        count = count_fasta(file)
        lines.append([file, 'fasta', count])
    for file in args.variant_read_ids:
        count = count_variant_read_ids(file)
        lines.append([file, 'variant_tsv', count])
    for file in args.pivoted_tsv_files:
        count = count_pivoted_tsv(file)
        lines.append([file, 'pivoted_tsv', count])
    for file in args.distinct_read_counts_files:
        count = count_pivoted_tsv(file)
        lines.append([file, 'distinct_read_counts', count])

    # print if no output file specified
    if args.output is None:
        for line in lines:
            print('\t'.join(map(str, line)))
        return

    # write the lines to the output file
    with open(args.output, 'wt') as f:
        for line in lines:
            f.write('\t'.join(map(str, line)) + '\n')

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

def count_variant_read_ids(file):

    # iterate over lines of file
    with use_open(file, 'rt') as f:
        count = sum(1 for line in f if line != '')

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
