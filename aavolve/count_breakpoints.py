#!/usr/bin/env python3

# count breakpoints in shuffled capsid reads using the assigned parents
# from assign_parents.py

import argparse
import csv

from aavolve.utils import use_open, sort_var_names

def main():
    
    parser = argparse.ArgumentParser('Count breakpoints in shuffled capsid reads')
    parser.add_argument('-i', '--input', help='Input file, most likely parents from assign_parents.py', required=True)
    parser.add_argument('-o', '--output', help='Output file, locations of breakpoints in each read', required=True)
    parser.add_argument('-s1', '--summary1', help='Output file, number of breakpoints in each read', required=True)
    parser.add_argument('-s2', '--summary2', help='Output file, number of breakpoints at each location', required=True)
    args = parser.parse_args()


    # count breakpoints
    n_break_reads, n_break_loc = count_breakpoints(args.input, args.output)

    # write summary files
    with use_open(args.summary1, 'wt') as out1:
        out1.write('read_id\tbreakpoints\n')
        for read, n_break in n_break_reads.items():
            out1.write(f'{read}\t{n_break}\n')

    with use_open(args.summary2, 'wt') as out2:
        out2.write('location\tbreakpoints\n')
        for loc, n_break in n_break_loc.items():
            out2.write(f'{loc}\t{n_break}\n')

def count_breakpoints(file, outfile):
    """
    For each line of the input file, identify breakpoints.
    Keep track of:
     - the number of breakpoints in each read
     - the number of breakpoints at each location
    """

    with use_open(file, 'rt') as handle, use_open(outfile, 'wt') as out:

        # open file using DictReader to get column names
        reader = csv.DictReader(handle, delimiter='\t')

        # sort column names
        variants = sort_var_names(reader.fieldnames[1:])
        
        # open output file using DictWriter
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()

        # for counting breakpoints
        n_break_reads = dict()
        n_break_loc = {col: 0 for col in variants}

        # for each line
        for line in reader:
            breaks = find_breakpoints(line, variants)

            # update breakpoint counts
            n_break_reads[line['read_id']] = len(breaks)

            for col in breaks:
                n_break_loc[col] += 1

            # write line to output file
            write_line = {'read_id': line['read_id']}
            for col in variants:
                write_line[col] = 1 if col in breaks else 0
            writer.writerow(write_line)

    return n_break_reads, n_break_loc

def find_breakpoints(line, sorted_cols):
    """
    Identify breakpoints in a line from file
    """

    breaks = []
    if len(sorted_cols) == 0:
        return breaks
    curr = get_var_set(line[sorted_cols[0]])
    
    # iterate over sorted columns
    for col in sorted_cols[1:]:

        line_vars = get_var_set(line[col])
        intersection = curr.intersection(line_vars)

        # if the current column is different from the previous
        if not intersection:

            # add the current column to the list of breakpoints
            breaks.append(col)

            # update the current column
            curr = line_vars
        
    return breaks


def get_var_set(vars):
    var_set = set(vars.split(','))
    assert 'NA' not in var_set
    return var_set


if __name__ == '__main__':
    main()
