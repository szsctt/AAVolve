#!/usr/bin/env python3

# Apply variants from a file to a reference
# If only variants in wide format are supplied, assume these are sequences to be applied
# If assigned parents (wide format) and parental variant sequences (long format) are supplied, get sequences 
# and apply to reference

import argparse
import csv
import tempfile
from Bio import Seq
from Bio import SeqIO

from scripts.utils import use_open, sort_var_names, get_parents, get_reference_name

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
    if args.parents:
        f = create_temp_seq_file(args.variants, args.parents)
        variants_file = f.name
    else:
        variants_file = args.variants
    
    # Open variants file
    with use_open(variants_file, 'rt', newline='') as f, use_open(args.output, 'wt', newline='') as o:
        reader = csv.DictReader(f, delimiter='\t')
        for i, row in enumerate(reader):
            _, seq = apply_variants(ref, row) 
            if args.translate:
                seq = Seq.translate(seq)
            
            if args.fasta:
                o.write(f">{i}\n{seq}\n")
            else:
                o.write(f"{seq}\n")

    if args.parents:
        f.close()

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
    offset = 0 #offset to account for prior indels
    last_pos = -1
    # iterate over variants
    for var in variants:
        pos, vartype = var.split(':')

        # substitution
        if vartype == 'sub':
            pos = int(pos)  + offset

            # don't apply variants that are before last_pos (e.g. deleted positions)
            if pos < last_pos:
                continue

            seq = seq[:pos] + row[var] + seq[pos+1:]
            last_pos = pos

        # insertion
        elif vartype == 'ins':

            # insertion with value '.' means no insertion
            if row[var] == '.':
                continue

            # get position
            pos = int(pos) + offset

            # don't apply variants that are before last_pos (e.g. deleted positions)
            if pos < last_pos:
                continue

            # apply insertion
            seq = seq[:pos] + row[var] + seq[pos:]
            offset += len(row[var])
            last_pos = pos

        # deletion`
        elif vartype == 'del':
            # get position
            start, end = pos.split("_")
            start = int(start) + offset
            end = int(end) + offset

            # deletion that doesn't have value '.' means no deletion
            if row[var] != '.':
                continue

            if start < last_pos:
                continue

            # extra bases at the end of the read (beyond the reference)
            # show up as insertions.  Exclude these
            if start >= len(ref):
                continue

            # apply deletion
            seq = seq[:start] + seq[end:]
            offset -= end - start
            last_pos = end
        
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

def create_temp_seq_file(variants, parents):
    """
    Create a temporary file with sequences for each variant
    """

    # read in parents
    vars = get_parents(variants)
    ref_name = get_reference_name(variants)

    f = tempfile.NamedTemporaryFile(mode='w+t')
    with use_open(parents, 'rt') as p:
        reader = csv.DictReader(p, delimiter='\t')

        # write header to output
        f.write('\t'.join(reader.fieldnames) + '\n')

        # iterate over rows
        for row in reader:
            
            # get count and read_id
            if 'count' in row:
                count = row['count']
                del row['count']
            if 'read_id' in row:
                read_id = row['read_id']
                del row['read_id']
            
            # iterate over remaining seuqnces
            seqs = []
            for pvar in row:

                if pvar not in vars:
                    raise ValueError(f"Variant {pvar} not found in parent variants file")
                
                # get possible parents for this variant
                parents = row[pvar].split(',')
                pos_vars = vars[pvar] # get the variant objects for this position

                # get the sequences of each parent allele
                par_seqs = []
                for i in parents:
                    # reference is not included in parents dict, so
                    # just get ref allele from any parent
                    if i == ref_name:
                        var = list(pos_vars.values())[0]
                        par_seqs.append(var.refbases())
                    # get alt allele
                    elif i in pos_vars:
                        par_seqs.append(pos_vars[i].qbases())
                    # otherwise, assume it's the reference allele
                    else:
                        pass

                assert len(par_seqs) > 0
                assert all([i==par_seqs[0] for i in par_seqs]) # all bases should be the same
                seqs.append(par_seqs[0])

            # write row to file
            if 'read_id' in reader.fieldnames:
                seqs.insert(0, read_id)
            if 'count' in reader.fieldnames:
                seqs.insert(0, count)
            f.write('\t'.join(seqs) + '\n')
  
    f.seek(0)
    return f





if __name__ == '__main__':
    main()
