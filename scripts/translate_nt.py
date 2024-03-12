import csv
import argparse

from Bio import Seq

from scripts.utils import use_open

def main():
    
    parser = argparse.ArgumentParser(description="Translate counts from nucleotide to amino acid")
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()

    with use_open(args.input, 'rt') as infile, use_open(args.output, 'wt') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=reader.fieldnames)

        # write header
        writer.writeheader()

        # translate 'sequence' column in each row
        for row in reader:
            assert 'sequence' in row 
            row['sequence'] = Seq.translate(row['sequence'])
            writer.writerow(row)

if __name__   == "__main__":
    main()
