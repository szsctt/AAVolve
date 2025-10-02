#!/usr/bin/env python3
"""
Plot frequency polynomials of repeat counts for tidehunter vs C3POa outputs.

Usage example:
python3 tidehunter_test/scripts/plot_repeats_comparison.py \
  --tidehunter tidehunter_test/output/tidehunter/tidehunter_output.fa \
  --c3poa tidehunter_test/output/c3poa/c3poa_output.fa \
  --out tidehunter_test/output/repeats_compare.png \
  --degree 6 --bins 40
"""
from __future__ import annotations

import argparse
import re
from collections import Counter
from typing import Iterable, List, Tuple

import numpy as np
import matplotlib.pyplot as plt

# Make sure the repository root is on sys.path so `aavolve` is importable when
# running this script directly (e.g. python3 tidehunter_test/scripts/plot_repeats.py ...)
import sys
from pathlib import Path

_this_file = Path(__file__).resolve()
_repo_root = None
for p in [_this_file] + list(_this_file.parents):
    if (p / 'aavolve').is_dir():
        _repo_root = p
        break
if _repo_root is None:
    # fallback to two levels up (typical layout)
    _repo_root = _this_file.parents[2]
sys.path.insert(0, str(_repo_root))

from aavolve.utils import get_repeats_from_r2c2_name

# --- helpers ----------------------------------------------------------------- #
def fasta_reader(path: str) -> Iterable[Tuple[str, str]]:
    """Simple FASTA reader: yield (header, seq) with seq as a single string.
    Raises on malformed FASTA.
    """
    seq_lines: List[str] = []
    header: str | None = None
    open_fn = open
    if path.endswith(".gz"):
        import gzip
        open_fn = gzip.open

    with open_fn(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                if header is None:
                    raise ValueError(f"FASTA file {path} contains sequence line before header")
                seq_lines.append(line.strip())
        if header is not None:
            yield header, "".join(seq_lines)


def parse_tidehunter_repeats(header: str) -> float:
    """Extract repeat count from tidehunter header.
    Expected pattern (example):
      >SRR29543177.5_rep0_3.0 7073_78_7020_...
    This extracts the numeric token after '_rep<id>_'.
    Raises ValueError if not found.
    """
    m = re.search(r"_rep\d+_([0-9]+(?:\.[0-9]+)?)", header)
    if not m:
        raise ValueError(f"Cannot parse repeats from tidehunter header: '{header}'")
    return float(m.group(1))


def extract_repeats_from_fasta(path: str, source: str) -> List[float]:
    """Extract numeric repeat counts from a FASTA for a given source type:
       source == 'tidehunter' -> use parse_tidehunter_repeats()
       source == 'c3poa' -> use aavolve.utils.get_repeats_from_r2c2_name()
    Raises on any parse error.
    """
    repeats = []
    for header, _seq in fasta_reader(path):
        if source == "tidehunter":
            r = parse_tidehunter_repeats(header)
        elif source == "c3poa":
            # get_repeats_from_r2c2_name returns numeric value or raises on bad name
            # R2C2/C3POa uses '0' for single-pass reads; tidehunter uses '1' for single-pass.
            # Add 1 to R2C2 counts to make conventions comparable.
            r = get_repeats_from_r2c2_name(header)
            # ensure numeric (function may return int/float) and shift by +1
            r = float(r) + 1.0
        else:
            raise ValueError("source must be 'tidehunter' or 'c3poa'")
        repeats.append(r)
    if not repeats:
        raise ValueError(f"No repeats extracted from {path} ({source})")
    return repeats


def frequency_and_polygon(repeats: List[float], bins: int | None):
    """Compute histogram frequencies and return centers and counts for a frequency polygon.
    Returns (centers, counts)
    """
    # Use integer-centered bins where possible: compute histogram over range
    min_r = int(np.floor(min(repeats)))
    max_r = int(np.ceil(max(repeats)))
    # create bin edges: if bins is None, create integer-centered bins for each integer repeat
    if bins is None:
        # one bin per integer repeat value from min_r..max_r inclusive
        bin_edges = np.arange(min_r - 0.5, max_r + 0.5 + 1.0, 1.0)
        if bin_edges.size < 2:
            bin_edges = np.array([min_r - 0.5, max_r + 0.5])
    else:
        bin_edges = np.linspace(min_r - 0.5, max_r + 0.5, max(2, bins))
    counts, edges = np.histogram(repeats, bins=bin_edges)
    # centers
    centers = (edges[:-1] + edges[1:]) / 2.0
    return centers, counts


# --- main -------------------------------------------------------------------- #
def get_args(argv=None):
    p = argparse.ArgumentParser(description="Compare repeat-count distributions for tidehunter and C3POa")
    p.add_argument("--tidehunter", required=True, help="tidehunter FASTA (two-line per entry) path")
    p.add_argument("--tidehunter-noargs", required=False, help="Optional tidehunter FASTA produced without -5/-3 args")
    p.add_argument("--c3poa", required=True, help="C3POa FASTA path")
    p.add_argument("--out", required=True, help="Output PNG path")
    p.add_argument("--bins", type=int, default=None, help="Number of bins used for histogram (default: one bin per integer repeat)")
    p.add_argument("--log-x", action='store_true', help="Use log scale on x-axis (requires positive bin centers)")
    return p.parse_args(argv)


def main(argv=None):
    args = get_args(argv)

    # extract repeats
    tide_repeats = extract_repeats_from_fasta(args.tidehunter, "tidehunter")
    tide_noargs_repeats = None
    if getattr(args, 'tidehunter_noargs', None):
        tide_noargs_repeats = extract_repeats_from_fasta(args.tidehunter_noargs, "tidehunter")
    c3poa_repeats = extract_repeats_from_fasta(args.c3poa, "c3poa")

    # frequency polygon data
    t_centers, t_counts = frequency_and_polygon(tide_repeats, bins=args.bins)
    c_centers, c_counts = frequency_and_polygon(c3poa_repeats, bins=args.bins)
    if tide_noargs_repeats is not None:
        tn_centers, tn_counts = frequency_and_polygon(tide_noargs_repeats, bins=args.bins)
    else:
        tn_centers, tn_counts = None, None

    # plot
    plt.figure(figsize=(9, 5))
    # draw frequency polygon (lines through bin centers)
    # vertical reference at 3 repeats
    try:
        plt.axvline(x=3, color='k', linestyle=':', linewidth=1)
    except Exception:
        pass
    if args.log_x:
        # ensure centers are positive; if zeros or negatives exist, shift all centers by +1
        all_centers = [c for c in [t_centers, c_centers, tn_centers] if c is not None]
        min_center = min(min(x) for x in all_centers)
        shift = 0.0
        if min_center <= 0:
            shift = 1.0 - min_center
            t_centers = t_centers + shift
            c_centers = c_centers + shift
            if tn_centers is not None:
                tn_centers = tn_centers + shift
        plt.xscale('log')
        plt.plot(t_centers, t_counts, color="C0", lw=2, marker='o', label="tidehunter freq poly", alpha=0.8)
        if tn_centers is not None:
            plt.plot(tn_centers, tn_counts, color="C2", lw=2, marker='o', label="tidehunter_noargs freq poly", alpha=0.7)
        plt.plot(c_centers, c_counts, color="C1", lw=2, marker='o', label="C3POa freq poly", alpha=0.7)
        # annotate if we shifted
        if shift > 0:
            plt.xlabel(f"Repeat count (shifted by +{shift:.1f} for log scale)")
        else:
            plt.xlabel("Repeat count")
    else:
        plt.plot(t_centers, t_counts, color="C0", lw=2, marker='o', label="tidehunter freq poly", alpha=0.8)
        if tn_centers is not None:
            plt.plot(tn_centers, tn_counts, color="C2", lw=2, marker='o', label="tidehunter_noargs freq poly", alpha=0.7)
        plt.plot(c_centers, c_counts, color="C1", lw=2, marker='o', label="C3POa freq poly", alpha=0.7)
    # ylabel/title/legend
    plt.ylabel("Frequency (reads)")
    plt.title("Repeat-count distribution: tidehunter vs C3POa")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    plt.close()


if __name__ == "__main__":
    main()