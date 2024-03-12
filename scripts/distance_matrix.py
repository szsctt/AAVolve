#!/usr/bin/env python3

# make a distnace matrix from amino acid sequences
# use Levenshtein distance to calculate the distance between two sequences

import argparse
import tempfile
import random
import csv
import subprocess
import io

import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator

import seaborn as sns
import matplotlib.pyplot as plt

from scripts.utils import use_open

def main():

    parser = argparse.ArgumentParser(description='make a distance matrix from amino acid sequences')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-p', '--plot', help='plot file')
    parser.add_argument('-d', '--distance-metric', help='distance metric', default='identity', choices=['identity', 'blosum62'])
    parser.add_argument('-s', '--max-seqs', help="Maximum number of sequences to read", type=int, default=1000)
    parser.add_argument('-x', '--selection', help='selection', default='random', choices=['random', 'first', 'last'])
    args = parser.parse_args()

    # read input file
    seqs = read_input(args.input, args.max_seqs, args.selection)

    if len(seqs) == 0:
        print('No sequences to process')
        # touch output files
        open(args.output, 'w').close()
        open(args.plot, 'w').close()
        return

    # calculate distance matrix
    mat = dmat(seqs, args.distance_metric)

    # write distance matrix
    np.savetxt(args.output, mat, delimiter='\t', fmt='%.10f') 

    # make and save plot
    if args.plot:
        sns.heatmap(mat, cmap='viridis')
        plt.savefig(args.plot)
    
def alignment(seqs):

    with tempfile.NamedTemporaryFile(mode='w') as handle:

        if len(seqs) == 0:
            return MultipleSeqAlignment([])

        # write sequences to a temporary fasta file
        for i, seq in enumerate(seqs):
            handle.write(f'>seq_{i}\n{seq}\n')

        handle.flush()

        # run mafft
        mafft_cline = ["mafft", '--auto', handle.name]
        r = subprocess.run(mafft_cline, capture_output=True)
        try:
            r.check_returncode()
        except subprocess.CalledProcessError as e:
            print(e.stderr.decode())
            raise e
        stdout = r.stdout.decode()

        # read with AlignIO
        return AlignIO.read(io.StringIO(stdout), 'fasta')


def dmat(seqs, metric):

    # do multiple sequence alignment
    aln = alignment(seqs)
    
    # calculate distance matrix
    calculator = DistanceCalculator(metric)
    dm = calculator.get_distance(aln)

    # convert to numpy array as lower triangular matrix
    dist = np.array([i + [0]*(len(aln) - len(i)) for i in dm.matrix])
    
    # convert lower triangular to full matrix
    return dist + dist.T


def read_input(infile, max_seqs, selection):

    # check inputs
    if max_seqs < 0:
        raise ValueError('max_seqs must be non-negative')
    if selection not in ['random', 'first', 'last']:
        raise ValueError("selection must be 'random', 'first', or 'last'")

    # read file
    with use_open(infile, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        seqs = []

        # skip header
        try:
            next(reader)
        except StopIteration:
            raise ValueError('Empty input file')
        
        # read rows
        for row in reader:
            if len(row) > 1:
                seqs.append(row[-1])
            else:
                seqs.append(row[0])

    print(f'Read {len(seqs)} sequences')
    if len(seqs) <= max_seqs:
        return seqs

    # downsample using specified method
    if selection == 'random':
        seqs = random.sample(seqs, max_seqs)
        print(f'Randomly selected {max_seqs} sequences')
    elif selection == 'first':
        seqs = seqs[:max_seqs]
        print(f'Selected first {max_seqs} sequences')
    elif selection == 'last':
        seqs = seqs[-max_seqs:]
        print(f'Selected last {max_seqs} sequences')

    return seqs


if __name__ == "__main__":
    main()