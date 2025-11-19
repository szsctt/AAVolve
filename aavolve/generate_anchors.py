#!/usr/bin/env python3

"""Generate 5' and 3' anchor sequences for a sample."""

import argparse
import random
import sys
from pathlib import Path

from aavolve.utils import use_open

NUCLEOTIDES = ("A", "C", "G", "T")


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generate 5' and 3' anchors of specified length")
    parser.add_argument("--length", "-l", type=int, required=True, help="Length of each anchor (5' and 3')")
    parser.add_argument("--output", "-o", required=True, help="Output FASTA file to write anchors to")
    parser.add_argument("--seed", help="Optional seed for reproducible anchor generation")
    return parser.parse_args(argv)


def write_anchors(length, output_path, seed=None):
    if length <= 0:
        raise ValueError("Anchor length must be a positive integer")

    output_path = Path(output_path)
    rng = random.Random(seed)
    five_prime = "".join(rng.choices(NUCLEOTIDES, k=length))
    three_prime = "".join(rng.choices(NUCLEOTIDES, k=length))

    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    with use_open(str(output_path), "wt") as handle:
        handle.write(
            f">5prime_anchor\n{five_prime}\n>3prime_anchor\n{three_prime}\n"
        )


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    write_anchors(args.length, args.output, seed=args.seed)


if __name__ == "__main__":
    main()
