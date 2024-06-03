#!/usr/bin/env python3

# comvert a table contaning each variant in an individual row
# to a table where each row is a read, with each variant in a column
# columns with no variant are filled with the reference allele
# and variants that don't match any parents are filled with "NA"

import argparse
import csv
from sys import argv

from aavolve.utils import get_variant, use_open, sort_var_names, get_reference_name, get_parents, make_var_groups

def main(sys_argv):

    args = get_args(sys_argv)

    # read in the parents file
    parents = get_parents(args.parents)

    # pivot the reads
    pivot_reads(args.input, args.read_ids, args.output_parents, args.output_seq, parents, args.remove_na, args.group_vars, args.group_dist, args.max_distance_frac)

def get_args(sys_argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help="input file", required=True)
    parser.add_argument('--read-ids', '-r', help='read ids file', required=True)
    parser.add_argument("--parents", "-p", help="parents file", required=True)
    parser.add_argument("--remove-na", action="store_true", help = "Remove variants that don't match any parents")
    parser.add_argument("--output-parents", "-o", help="output file with possible parents for each variant", required=True)
    parser.add_argument("--output-seq", "-O", help="output file with sequence for each variant", required=True)
    parser.add_argument("--group-vars", '-g', action="store_true", help="Group adjacent variants and use the parent with the lowest hamming distance")
    parser.add_argument("--group-dist", help="Variants separated by at most this number of nucleotides will be grouped for assigning parents, and lowest hamming distance parent will be used", default=1, type=int)
    parser.add_argument("--max-distance-frac", help="When grouping variants, if distance divided by variant number in group is less than this for a parent, the parent will be assigned to the group", type=float, default=0.2)
    args = parser.parse_args(sys_argv)
    return args

def pivot_reads(infile, in_read_ids, outfile_parents, outfile_seq, parents, remove_na, group, group_dist, max_dist_frac):
    """
    Write one line per read to outfile
    Columns in outfile are read_id, variant1, variant2, variant3, ...
    """

    # get distinct variant ids from parents
    parent_var_ids = list(set(parents.keys()))

    # sort variant ids by position
    parent_var_ids = sort_var_names(parent_var_ids)

    # create a header for writing later
    header = ["read_id"] + list(parent_var_ids)

    # group parent var ids
    parent_var_ids = make_var_groups(parent_var_ids, group, group_dist)

    # get names of parents, including reference
    wt_name = get_reference_name(infile)
    parent_names = set([wt_name])
    for d in parents.values():
        for k in d.keys():
            parent_names.add(k)

    # rearrange parents to get dict with variants from each parent as values
    all_parent_group_vars = {}
    # iterate over groups
    for group in parent_var_ids:
        all_parent_group_vars[group] = {}
        # for each group, get variants
        for parent_name in parent_names:
            all_parent_group_vars[group][parent_name] = {parents[var][parent_name].var_id():parents[var][parent_name] for var in group if parent_name in parents[var]}
    assert len(all_parent_group_vars) == len(parent_var_ids)

    # open result files for writing
    with use_open(outfile_parents, "wt") as par, use_open(outfile_seq, "wt") as seq:
        # write header for parents file
        writer_par = csv.writer(par, delimiter="\t")
        writer_par.writerow(header)
            
        # write header for seq file
        writer_seq = csv.writer(seq, delimiter="\t")
        writer_seq.writerow(header)

        # iterate over reads
        for read_id, read in get_reads(infile, in_read_ids):
            
            # start new row with read id
            row_seq, row_par = [read_id], [read_id]
            
            # iterate over possible parental variants
            for group in parent_var_ids:
                
                # collect variants for this group
                group_vars = [read[var] for var in group if var in read.keys()]

                # find closest parent from parent_group_vars -> returns dict with parent as key and variants as value
                # there may be more than one parent with the same distance
                group_pars = closest_parent(group, group_vars, all_parent_group_vars[group], max_dist_frac)

                # if we have multiple closest parents with different variants, pick one randomly for the sequence file
                group_par_names = ','.join(sorted(list(group_pars.keys())))
                group_par = set(group_pars.keys()).pop()
                
                # if parent is NA, no closest parent could be identified, so just add NA to rows
                if group_par == 'NA':
                    for var in group:
                        row_seq.append('NA')
                        row_par.append('NA')
                    continue
                
                # get variants for parent
                parent_group_vars = group_pars[group_par]

                # assign variants and parent for each variant in group
                for var in group:

                    # add parent name(s)
                    row_par.append(group_par_names)
                    
                    # if variant not in parent, add wild-type sequence
                    if var not in parent_group_vars.keys():
                        # get reference sequence for this variant from another parent
                        bases = ''
                        for parent in parent_names:
                            if var in all_parent_group_vars[group][parent]:
                                bases = all_parent_group_vars[group][parent][var].refbases()
                                break
                        assert bases != ''
                        row_seq.append(bases)
                    # add variant sequence
                    else:
                         row_seq.append(parent_group_vars[var].qbases())

            # check if any variants didn't match any parents
            if remove_na and "NA" in row_par:
                continue

            # write row to outfile
            assert len(header) == len(row_par)
            assert len(row_par) == len(row_seq)
            writer_par.writerow(row_par)
            writer_seq.writerow(row_seq)

def closest_parent(group, read_vars, parents, max_dist_frac):

    dists = {}
    # get distance for each parent
    for par in parents.keys():
        # get all the variants in the read and the parent
        read_variants = set(i.var_id() for i in read_vars)
        parent_variants = set(i for i in parents[par].keys())
        dist = 0
        # for each variant, distance is incremented if the variant is different between the read and parent
        for var in group:
            # both read and parent have this variant
            if var in read_variants and var in parent_variants:
                # but they're different
                if parents[par][var] not in read_vars:
                    dist += 1 
            else:
                # variant is in read but not parent or vice versa
                if var in read_variants or var in parent_variants:
                    dist += 1
        dists[par] = dist

    # check if any distances are below maximum
    max_dist = int(len(group)* max_dist_frac)
    if any(i <= max_dist for i in dists.values()):
        par = min(dists, key=dists.get)
        par = [i for i in dists.keys() if dists[i] == dists[par]] # get any tied parents
    else:
        return {'NA': 'NA'} # if all distances are above max, return NAs

    return {i: parents[i] for i in par}

def collect_read_vars(file_name):

    with use_open(file_name, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        # collect variants for one read at a time
        rid, vars = '', {}
        for row in reader:
            
            # if this is a new read, yield
            if row['query_name'] != rid:
                if len(vars) > 0:
                    yield rid, vars
                rid, vars = row['query_name'], {}

            # get variants from line
            var = get_variant(row)
            vars[var.var_id()] = var

    # yield the last read
    if len(vars) > 0:
        yield rid, vars

def get_read_id(file_name):

    with open(file_name, "rt") as handle:
        for line in handle:
            rid = line.strip()
            if rid != '':
                yield rid

def get_reads(file_name, in_read_ids):
    """
    Collect the variants from each read in the input file
    """

    reader = collect_read_vars(file_name)
    read_ids = get_read_id(in_read_ids)
    var_read_id, vars, rid = "", {}, ""

    for var_read_id, vars in reader:
        # get read id from read id file
        rid = next(read_ids)

        # check to see if this read id is in the input file
        while var_read_id != rid:
            yield rid, {}
            rid = next(read_ids)
        yield var_read_id, vars
    
    # yield any read ids after end of reads with variants
    try:
        rid = next(read_ids)
    except StopIteration:
        return

    while var_read_id != rid:
        yield rid, {}
        try:
            rid = next(read_ids)
        except StopIteration:
            return

if __name__ == "__main__":
    main(argv[1:])

