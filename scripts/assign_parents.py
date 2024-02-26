#!/usr/bin/env python3

# assign the most likely parents for each variant
# by contintuing through the read until the next variant doesn't match
# any of the current parents

import argparse
import csv

from scripts.utils import use_open

def main():
    
    parser = argparse.ArgumentParser(description='assign the most likely parents for each variant')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()

    # open files
    with use_open(args.input, 'rt') as infile, use_open(args.output, 'wt') as outfile:
       
        # create csv reader and writer
        reader = csv.reader(infile, delimiter='\t') 
        writer = csv.writer(outfile, delimiter='\t')

        # get header
        header = next(reader)

        # write header
        writer.writerow(header)

        for row in reader:
            # get parents for each variant
            row = get_parents(row)
            # write row
            writer.writerow(row)


def get_parents(row):
    # get the parents for each variant
    # by contintuing through the read until the next variant doesn't match
    # any of the current parents

    # return if no variants
    if len(row) == 1:
        return row
    
    assigned = row[1:].copy()
    
    curr_parents = set(assigned[0].split(","))
    last_breakpoint = 0

    for i, entry in enumerate(assigned):

        # get parents for this entry
        entry_parents = set(entry.split(","))

        # if NA, skip
        if entry_parents == {"NA"}:
            continue

        # get there is an intersection between current parents and entry parents
        intersection = curr_parents.intersection(entry_parents)
        if intersection:
            # if there is an intersection, update the current parents
            curr_parents = intersection
        else:
            # for writing, sort parents in alphabetical order
            write_parents = ','.join(sorted(list(curr_parents)))

            # assign the current parants to all positions between here
            # and last breakpoint
            assigned[last_breakpoint:i] = [write_parents for _ in range(last_breakpoint, i)]

            # reset current parents
            curr_parents = entry_parents

            # assign last breakpoint
            last_breakpoint = i

    # assign last lot of parents
    write_parents = ','.join(sorted(list(curr_parents)))
    assigned[last_breakpoint:] = [write_parents for _ in range(last_breakpoint, len(assigned))]

    return row[0:1] + assigned
    
if  __name__ == '__main__':
    main()
