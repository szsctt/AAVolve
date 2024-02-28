#!/usr/bin/env python3

# Apply variants from a file to a reference
# If only variants in wide format are supplied, assume these are sequences to be applied
# If assigned parents (wide format) and parental variant sequences (long format) are supplied, get sequences 
# and apply to reference

import argparse
import csv
from Bio import Seq
from Bio import SeqIO

from scripts.utils import use_open, sort_var_names

def main():
    
    parser = argparse.ArgumentParser(description='Apply variants (wide format) to a reference sequences (fasta) to generate a table containing count and sequence')
    parser.add_argument('-v', '--variants', help='Variants file (wide format)', required=True)
    parser.add_argument('-p', '--parents', help='Parental variants (long format)')
    parser.add_argument('-r', '--reference', help='Reference sequence (fasta)', required=True)
    parser.add_argument('-t', '--translate', help='Translate to amino acids', action='store_true')
    parser.add_argument('-f', '--fasta', help='Output in fasta format instead of tab-delimited', action='store_true')
    parser.add_argument('-o', '--output', help='Output file', required=True)
    args = parser.parse_args()

    # Read reference sequence
    ref = read_reference(args.reference)

    # if parents are supplied, read in and use to create temp file that contains sequences for each variant
    # TODO
    
    # Open variants file
    with use_open(args.variants, 'rt', newline='') as f, use_open(args.output, 'wt', newline='') as o:
        reader = csv.DictReader(f, delimiter='\t')
        for i, row in enumerate(reader):
            _, seq = apply_variants(ref, row) 
            if args.translate:
                seq = Seq.translate(seq)
            
            if args.fasta:
                o.write(f">{i}\n{seq}\n")
            else:
                o.write(f"{seq}\n")

def apply_variants(ref, row):
    """
    Apply variant to reference sequence
    """

    # Count is 1 if not specified
    if 'count' in row:
        count = row['count']
        del row['count']
    else:
        count = 1
        if 'read_id' in row:
            del row['read_id']

    # order variants by position
    variants = sort_var_names(row.keys())
    
    seq = ref
    # iterate over variants
    for var in variants:
        pos, vartype = var.split(':')

        # substitution
        if vartype == 'sub':
            pos = int(pos)  
            seq = seq[:pos] + row[var] + seq[pos+1:]

        # insertion
        elif vartype == 'ins':

            # insertion with value '.' means no insertion
            if row[var] == '.':
                continue

            # apply insertion
            pos = int(pos)
            seq = seq[:pos] + row[var] + seq[pos:]

        # deletion`
        elif vartype == 'del':
            # get position
            start, end = pos.split("_")
            start = int(start)
            end = int(end)

            # deletion that doesn't have value '.' means no deletion
            if row[var] != '.':
                continue

            # extra bases at the end of the read (beyond the reference)
            # show up as insertions.  Exclude these
            if start >= len(ref):
                continue

            # apply deletion
            seq = seq[:start] + seq[end:]
        
        else:
            raise ValueError(f"Unknown variant type: {vartype}")      

    return count, seq


def read_reference(reference):
    """
    Read reference sequence from fasta file
    """
    recs = [i for i in SeqIO.parse(reference, "fasta")]
    assert len(recs) == 1
    return str(recs[0].seq)

if __name__ == '__main__':
    main()
