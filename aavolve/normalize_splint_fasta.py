#!/usr/bin/env python3

import argparse
from sys import argv

from aavolve.utils import seq_generator, use_open


def normalize_splint_fasta(in_file: str, out_file: str) -> None:
    records = []
    with use_open(in_file, "rt") as handle:
        for name, seq in seq_generator(handle):
            records.append((name, seq))

    if len(records) != 1:
        raise ValueError(
            f"Splint FASTA must contain exactly 1 sequence, found {len(records)} in {in_file}"
        )

    (_, seq) = records[0]

    with open(out_file, "w") as out:
        out.write(">splint\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i : i + 60] + "\n")


def get_args(args):
    parser = argparse.ArgumentParser(
        description="Ensure splint FASTA has one record and rename it to 'splint'."
    )
    parser.add_argument("--input", "-i", required=True, help="Input splint FASTA (optionally gzipped)")
    parser.add_argument("--output", "-o", required=True, help="Output normalized splint FASTA")
    return parser.parse_args(args)


def main(args):
    arg = get_args(args)
    normalize_splint_fasta(arg.input, arg.output)


if __name__ == "__main__":
    main(argv[1:])
