#!/usr/bin/env python3

"""Apply stored anchor sequences to FASTA/FASTQ inputs."""

import argparse
import sys
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq

from pathlib import Path

from aavolve.utils import use_open

DEFAULT_PHRED = 40
ANCHOR_IDS = {"5prime_anchor", "3prime_anchor"}


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Add 5' and 3' anchors to sequences")
    parser.add_argument("--input", "-i", required=True, help="Input FASTA/FASTQ file")
    parser.add_argument("--output", "-o", required=True, help="Output file with anchors applied")
    parser.add_argument("--anchors", "-a", required=True, help="FASTA file containing 5' and 3' anchors")
    return parser.parse_args(argv)


def load_anchors(path: str) -> Tuple[str, str]:
    anchors = {}
    with use_open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            anchors[record.id] = str(record.seq).upper()
    missing = ANCHOR_IDS - anchors.keys()
    if missing:
        raise ValueError(f"Anchor file {path} is missing entries for {', '.join(sorted(missing))}")
    return anchors["5prime_anchor"], anchors["3prime_anchor"]


def detect_format(path: str) -> str:
    with use_open(path, "rt") as handle:
        first_char = handle.read(1)
    if first_char == '>':
        return "fasta"
    if first_char == '@':
        return "fastq"
    raise ValueError(f"Unable to detect sequence format for {path!r}")


def pad_quality(record, qualities, front_len, back_len):
    anchor_quality = [DEFAULT_PHRED]
    record.letter_annotations["phred_quality"] = (
        anchor_quality * front_len + qualities + anchor_quality * back_len
    )


def apply_anchors(input_path: str, output_path: str, anchors_path: str):
    five_prime, three_prime = load_anchors(anchors_path)
    seq_format = detect_format(input_path)

    front_len = len(five_prime)
    back_len = len(three_prime)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with use_open(input_path, "rt") as in_handle, use_open(output_path, "wt") as out_handle:
        for record in SeqIO.parse(in_handle, seq_format):
            qualities = None
            if seq_format == "fastq":
                if "phred_quality" not in record.letter_annotations:
                    raise ValueError(f"FASTQ record {record.id} lacks phred_quality annotations")
                qualities = list(record.letter_annotations["phred_quality"])
            record.letter_annotations = {}
            combined = five_prime + str(record.seq).upper() + three_prime
            record.seq = Seq(combined)
            if seq_format == "fastq":
                pad_quality(record, qualities, front_len, back_len)
            SeqIO.write(record, out_handle, seq_format)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    apply_anchors(args.input, args.output, args.anchors)


if __name__ == "__main__":
    main()
