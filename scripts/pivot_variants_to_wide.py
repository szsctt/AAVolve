#!/usr/bin/env python3

# comvert a table contaning each variant in an individual row
# to a table where each row is a read, with each variant in a column
# columns with no variant are filled with the reference allele
# and variants that don't match any parents are filled with "NA"

import argparse
import csv
from sys import argv

from scripts.utils import get_variant, use_open, sort_var_names, get_reference_name, get_parents

def main(sys_argv):

    args = get_args(sys_argv)

    # read in the parents file
    parents = get_parents(args.parents)

    # pivot the reads
    pivot_reads(args.input, args.output_parents, args.output_seq, parents, args.remove_na)

def get_args(sys_argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help="input file", required=True)
    parser.add_argument("--parents", "-p", help="parents file", required=True)
    parser.add_argument("--remove-na", action="store_true", help = "Remove variants that don't match any parents")
    parser.add_argument("--output-parents", "-o", help="output file with possible parents for each variant", required=True)
    parser.add_argument("--output-seq", "-O", help="output file with sequence for each variant", required=True)
    args = parser.parse_args(sys_argv)
    return args

def pivot_reads(infile, outfile_parents, outfile_seq, parents, remove_na):
    """
    Write one line per read to outfile
    Columns in outfile are read_id, variant1, variant2, variant3, ...
    """

    # get distinct variant ids from parents
    parent_var_ids = list(set(parents.keys()))

    # sort variant ids by position
    parent_var_ids = sort_var_names(parent_var_ids)

    # get names of parents, including reference
    wt_name = get_reference_name(infile)
    parent_names = [wt_name]
    for d in parents.values():
        for k in d.keys():
            if k not in parent_names:
                parent_names.append(k)

    # write a header to outfile
    header = ["read_id"] + list(parent_var_ids)

    with use_open(outfile_parents, "wt") as par, use_open(outfile_seq, "wt") as seq:
        # write header for parents file
        writer_par = csv.writer(par, delimiter="\t")
        writer_par.writerow(header)
            
        # write header for seq file
        writer_seq = csv.writer(seq, delimiter="\t")
        writer_seq.writerow(header)

        # iterate over reads
        for read_id, read in get_reads(infile):
            
            # start new row with read id
            row_seq, row_par = [read_id], [read_id]
            
            # iterate over possible parental variants
            for parent_var_id in parent_var_ids:
                # check if this read has this variant
                if parent_var_id in read.keys():

                    # get variant from read
                    read_var = read[parent_var_id]

                    # get parents that have this variant
                    vars_parents = [k for k,v in parents[parent_var_id].items() if str(read_var) == str(v)]

                    # if there are any parents that have this variant
                    if len(vars_parents) > 0:
                        # add sequence to row
                        row_seq.append(read_var.qbases())
                        # add parent name to row
                        row_par.append(",".join(vars_parents))
        
                    # otherwise, add NA to row
                    else:
                        row_seq.append("NA") # could potentially be the alternate allele?
                        row_par.append("NA")

                # otherwise, sequence is wild-type
                else:

                    # get any variant from parents to get the reference allele
                    parent_var = list(parents[parent_var_id].values())[0]
                    # sequnece is reference
                    row_seq.append(parent_var.refbases())
                    # parents are all those that have the wild-type allele
                    alt_parents = set((*parents[parent_var_id].keys(), 'non_parental'))
                    vars_parents = [i for i in parent_names if i not in alt_parents]
                    row_par.append(",".join(vars_parents))

            # check if any variants didn't match any parents
            if remove_na and "NA" in row_par:
                continue

            # write row to outfile
            writer_par.writerow(row_par)
            writer_seq.writerow(row_seq)


def get_reads(file_name):
    """
    Collect the variants from each read in the input file
    """

    with use_open(file_name, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        
        # get first variant and read id
        try:
            line = next(reader)
        except StopIteration:
            return
        read_id = line['query_name']
        var = get_variant(line)
        buffer = {var.var_id(): var}
        
        for row in reader:

            # get variant from row
            var = get_variant(row)

            # check if we're still on the same read
            if row['query_name'] == read_id:
                buffer[var.var_id()] = var
            # if not, yield the read and start a new buffer
            else:
                yield read_id, buffer
                read_id = row['query_name']
                buffer = {var.var_id(): var}

        # yield the last read
        yield read_id, buffer

if __name__ == "__main__":
    main(argv[1:])

