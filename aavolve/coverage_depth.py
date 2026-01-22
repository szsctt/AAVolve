#!/usr/bin/env python3

import argparse
import os
import sys
from typing import Optional

import pysam

from aavolve.utils import use_open


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute per-position read depth across the reference.")
    parser.add_argument("--bam", required=True, help="Input BAM (indexed).")
    parser.add_argument("--reference-name", default=None, help="Reference/contig name (default: first BAM contig).")
    parser.add_argument("--output", "-o", required=True, help="Output TSV(.gz) with columns: ref, pos, depth.")
    return parser.parse_args(argv)


def _is_primary(read: pysam.AlignedSegment) -> bool:
    return read.is_mapped and not read.is_secondary and not read.is_supplementary


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        if args.reference_name is None:
            if not bam.references:
                raise ValueError(f"No references found in BAM {args.bam!r}")
            ref_name = bam.references[0]
        else:
            ref_name = args.reference_name
            if ref_name not in bam.references:
                raise ValueError(f"Reference {ref_name!r} not found in BAM {args.bam!r}")

        ref_len = bam.get_reference_length(ref_name)
        depths = [0] * ref_len

        for pileup_col in bam.pileup(
            ref_name,
            start=0,
            end=ref_len,
            truncate=True,
            stepper="nofilter",
            ignore_overlaps=False,
            ignore_orphans=False,
            min_base_quality=0,
        ):
            pos = pileup_col.reference_pos
            if pos < 0 or pos >= ref_len:
                continue
            depth = 0
            for pr in pileup_col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                read = pr.alignment
                if not _is_primary(read):
                    continue
                depth += 1
            depths[pos] = depth

    with use_open(args.output, "wt", newline="") as handle:
        handle.write("ref\tpos\tdepth\n")
        for i, d in enumerate(depths, start=1):
            handle.write(f"{ref_name}\t{i}\t{d}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

