#!/usr/bin/env python3

# the goal of this script is to calculate allele frequency.
# that is, the substitution T51G is distint from T51A

# always output high-frequency (above threshold args.frac) variants
# if parental variants are not provided, then output all high-frequency variants
# if parental variants are provided, then output high-frequency, non-parental variants

# optionally include the following outputs:
# 1. all frequencies
# 2. parental frequencies

# these frequencies can be used to extract high-frequency variants that can be combined
# with parental variants to account for high-frequency non-parental variants, and fed to
# filter_non_parental_variants.py to include these along with parental variants

# note that to to enable combination with the parental variants, I keep the same columns
# as this file.  This includes the reference and alternate allele 


import argparse

from scripts.utils import (get_variants_set, get_variant, get_reference_name,
                            use_open, read_variant_file, get_header)

def main():

    parser = argparse.ArgumentParser("Filter variants that didn't come from any of the parents")
    parser.add_argument("--input", "-i", help="Input variants", required=True)
    parser.add_argument("--parents", "-p", help="Parental variants to exclude (optional)")
    parser.add_argument("--output", "-o", help="Ouput file for high-frequency, non-parental variants", required=True)
    parser.add_argument("--output-all-freqs", "-oa", help="Output file for all frequencies")
    parser.add_argument("--output-parental-freqs", "-op", help="Output file for parental frequencies")
    parser.add_argument("--frac", "-f", help="Fraction of reads a variant must be seen in to be retained", type=float, default=0.1)
    args = parser.parse_args()
    
    # get parental variants
    if args.parents is not None:
        parents = get_variants_set(args.parents)
    else:
        parents = None
        
    # get frequency of each variant
    freqs = get_variant_frequency(args.input)
    
    # output frequencies of all variants
    if args.output_all_freqs is not None:
        if parents:
            write_freqs(freqs, parents, args.output_all_freqs, parents_only=False)
        else:
            write_freqs(freqs, {}, args.output_all_freqs, parents_only=False)
    
    # output frequencies of parental variants
    if args.output_parental_freqs is not None:
        if parents:
            write_freqs(freqs, parents, args.output_parental_freqs,parents_only=True) 
        else:
            write_freqs({}, {}, args.output_parental_freqs,parents_only=True)

    # filter by desired fraction
    freqs = {k:v for k, v in freqs.items() if v > args.frac}
    print(f'Found {len(freqs)} variants with a frequency above {args.frac}')
 
    # filter out parental variants
    if parents is not None:
        parents_str = set([str(p) for p in parents])
        freqs = {k:v for k, v in freqs.items() if not str(k) in parents_str}
        print(f'Found {len(freqs)} variants with a frequency above {args.frac} that are not parental variants')

    # output filtered variants
    write_variants(freqs, args.output, args.input)
    
def get_variant_frequency(file):
    """
    Calculate the frequency of each distinct variant
    Also keep track of reference and alternate alleles for each distinct variant
    """
    
    counts = dict()
    read_ids = set() # to count total number of reads
    
    # iterate over lines of file    
    for row in read_variant_file(file):
    
        # get read id
        read_ids.add(row['query_name'])
        
        # get variant id
        var = get_variant(row)
        
        # if we haven't seen this variant before
        if var not in counts.keys():
            counts[var] = 0
            
        # increment count
        counts[var] += 1
            
    # divide each count by library size to get fraction
    counts = {k: v/len(read_ids) for k, v in counts.items()}
         
    return counts

def write_freqs(freqs, parents,  filename, parents_only=False):
    """
    Write frequencies of parental variants to output file
    If parents_only is True, then only write frequencies of parental variants 

    """
    with use_open(filename, 'wt', newline='') as handle:
        
        # check for no frequencies
        if len(freqs) == 0:
            return
            
        # add header
        header = next(iter(freqs.keys())).header(shorter=True)
        header = header[:-1] + '\tfreq\n'
        handle.write(header)

        # write rows
        parents_str = set([str(p) for p in parents])
        for var, f in freqs.items():
                # is the variant parental or non-parental?
                par = "parental" if str(var) in parents_str else "non_parental"
                # skip if we only want parental variants
                if parents_only and par == "non_parental":
                    continue
                # write row
                row = var.print_line(par)[:-1]
                row = row + f"\t{f}\n"
                handle.write(row)



def write_variants(freqs, filename, input_file):
    """
    Write the high-frequency variants to file in same format as input file
    """  
    # need to match the format of parents_file in the output
    # so figure out if we used shorter or longer input

    header = get_header(input_file)
    if len(freqs) > 0:
        ref = get_reference_name(input_file)

    # open file and write header
    with use_open(filename, 'wt', newline='') as handle:
        handle.write(header)
        # iterate over variants
        for i, var in enumerate(freqs.keys()):
            # write line
            line = var.print_line(query_name=f'non_parental_{i}', ref_name = ref)
            handle.write(line)
             
         
    
if __name__=="__main__":
    main()
