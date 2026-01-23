#!/usr/bin/env python3

import argparse
import os
import sys
from typing import Optional

from aavolve.utils import read_variant_file


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Warn when high-frequency non-parental variants fall outside the parental-variant window "
            "while reads are not required to align end-to-end."
        )
    )
    parser.add_argument("--high-freq-variants", required=True, help="TSV/TSV.GZ of high-frequency non-parental variants.")
    parser.add_argument(
        "--parent-first-last",
        required=True,
        help="Two-line file containing the first and last parental variant positions (empty if one parent).",
    )
    parser.add_argument("--n-parents", required=True, type=int, help="Number of parent sequences in parent_file.")
    parser.add_argument("--include-non-parental", required=True, choices=("True", "False"))
    parser.add_argument("--require-end-to-end-alignment", required=True, choices=("True", "False"))
    parser.add_argument("--output", "-o", required=True, help="Output warning text file (empty if no warning).")
    return parser.parse_args(argv)


def _read_first_last(path: str) -> tuple[int, int]:
    with open(path, "rt") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"Expected two lines (first/last) in {path!r}, got {len(lines)}")
    return int(lines[0]), int(lines[1])


def _pos_range(pos_text: str) -> tuple[int, int]:
    if "_" in pos_text:
        start_text, end_text = pos_text.split("_", 1)
        return int(start_text), int(end_text)
    pos = int(pos_text)
    return pos, pos


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)

    include_non_parental = args.include_non_parental == "True"
    require_end_to_end_alignment = args.require_end_to_end_alignment == "True"

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    if not include_non_parental or require_end_to_end_alignment or args.n_parents <= 1:
        with open(args.output, "wt") as handle:
            handle.write("")
        return 0

    first, last = _read_first_last(args.parent_first_last)

    outside_positions: list[tuple[int, int]] = []
    for row in read_variant_file(args.high_freq_variants):
        pos_text = row.get("pos")
        if pos_text is None:
            continue
        start, end = _pos_range(pos_text)
        if end < first or start > last:
            outside_positions.append((start, end))

    if not outside_positions:
        with open(args.output, "wt") as handle:
            handle.write("")
        return 0

    min_outside = min(start for start, _ in outside_positions)
    max_outside = max(end for _, end in outside_positions)

    with open(args.output, "wt") as handle:
        handle.write(
            "High-frequency non-parental variants were detected outside the parental-variant window "
            f"[{first}, {last}] (outside range observed: [{min_outside}, {max_outside}], count={len(outside_positions)}).\n"
        )
        handle.write(
            "Reads were not required to align end-to-end (`require_end_to_end_alignment=False`), so allele frequencies "
            "near read ends may be diluted by partial-coverage reads.\n"
        )
        handle.write(
            "If end-of-read variants are important for your analysis, consider setting `require_end_to_end_alignment=True` "
            "(or trimming reads so alignments are end-to-end).\n"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

