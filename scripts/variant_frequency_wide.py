#!/usr/bin/env python3

# Calculate the frequency of variants for each column of file

import argparse
import csv

from scripts.utils import use_open

def main():
    parser = argparse.ArgumentParser(description='Calculate the frequency of variants for each column of file')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    parser.add_argument('--split-counts', help='Counts for mutiple parents are split between each parent (i.e. for AAV2,AAV8 counts for AAV2 and AAV8 are incremented by 0.5), default is not to split (so both incremented by 1)', action='store_true')
    args = parser.parse_args()


    freqs = count_freqs(args.input, args.split_counts)

    # write output
    with use_open(args.output, 'wt') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['variant', 'parent', 'frequency'])

        for var, counts in freqs.items():
            for parent, count in counts.items():
                writer.writerow([var, parent, count])

def count_freqs(input_file, split_counts=False):

    # get freqs from input file
    read_count = 0
    with use_open(input_file, 'rt') as in_file:
        reader = csv.DictReader(in_file, delimiter='\t')

        # check for empty file
        if reader.fieldnames is None:
            raise ValueError('No reads found in input file')
    
        # if no variants, return empty dict
        if len(reader.fieldnames) == 1:
            return {}

        # get number of fields
        num_fields = len(reader.fieldnames)
        
        # to store counts
        freqs = {variant: {} for variant in reader.fieldnames[1:]}

        # iterate through reads
        for row in reader:
            read_count += 1
            assert len(row) == num_fields, f"Row {row} has {len(row)} fields, expected {num_fields}"

            # for each variant, increment count for each parent
            for var, item in row.items():
                # skip read id colum
                if var == 'read_id':
                    continue
                
                # if only one parent, increment count
                if ','not in item:
                    if item not in freqs[var]:
                        freqs[var][item] = 0
                    freqs[var][item] += 1
                    continue

                # otherwise, split counts between parents
                parents = item.split(',')
                count_split = 1/len(parents)
                for parent in parents:
                    
                    # if parent not in freqs, add it
                    if parent not in freqs[var]:
                        freqs[var][parent] = 0
                    
                    # increment count
                    if split_counts:
                        freqs[var][parent] += count_split
                    else:
                        freqs[var][parent] += 1

    # calculate frequencies
    for var, counts in freqs.items():
        for parent, count in counts.items():
            freqs[var][parent] = count/read_count

    return freqs

if __name__ == "__main__":
    main()
