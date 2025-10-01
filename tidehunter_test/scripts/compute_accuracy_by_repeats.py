#!/usr/bin/env python3
"""Compute per-repeat-count accuracy for tidehunter and C3POa aligned reads.

Reads two BAM files (one per tool) aligned to the reference and computes
percent identity per read as (matches / aligned_bases) * 100. Groups reads by
their repeat count (extracted from the read header). For C3POa/R2C2 counts,
adds +1 to make them comparable to tidehunter convention.

Outputs a TSV with columns: tool, repeats, n_reads, mean_id, median_id, sd_id
Optionally writes a simple PNG comparing mean identity vs repeats.

Requires pysam and matplotlib (for plotting).
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from statistics import mean, median, pstdev
from typing import Dict, List, Tuple

try:
    import pysam
except Exception:
    pysam = None

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

# try to import the canonical parser from the package; provide a safe
# fallback when running scripts from the repo where the package may not be
# installed into the environment.
aav_utils = None
try:
    from aavolve import utils as aav_utils  # type: ignore
except Exception:
    try:
        # attempt to add the repo root to sys.path so local package imports work
        import os

        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)
        from aavolve import utils as aav_utils  # type: ignore
    except Exception:
        aav_utils = None


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Compute accuracy by repeat count from aligned BAMs")
    p.add_argument("--tide-bam", required=True, help="BAM file for tidehunter alignments")
    p.add_argument("--tide-bam-noargs", required=False, help="Optional BAM file for tidehunter run without -5/-3 arguments (compare as separate condition)")
    p.add_argument("--c3poa-bam", required=True, help="BAM file for C3POa/R2C2 alignments")
    p.add_argument("--out-tsv", required=True, help="Output TSV path")
    p.add_argument("--plot", help="Optional PNG path to plot median identity vs repeats")
    p.add_argument("--length-plot", help="Optional PNG path to plot median aligned length vs repeats (if omitted, derived from --plot by adding _length.png)")
    return p.parse_args(argv)


def get_repeats_from_tide_header(name: str) -> int:
    m = re.search(r"_rep\d+_([0-9]+(?:\.[0-9]+)?)", name)
    if not m:
        raise ValueError(f"Cannot parse repeat count from tidehunter header: {name}")
    return int(float(m.group(1)))


def get_repeats_from_c3poa_header(name: str) -> int:
    """Delegate R2C2/C3POa header parsing to aavolve.utils.get_repeats_from_r2c2_name.

    The helper implements the repository's canonical parsing and already handles
    any necessary adjustments (such as +1). Return an int repeat count or raise
    ValueError if parsing fails.
    """
    # Prefer the canonical helper when available
    hdr = name[1:] if name.startswith('>') else name
    val = aav_utils.get_repeats_from_r2c2_name(hdr)
    if val is None:
        raise ValueError(f"Cannot parse repeat count from c3poa header: {name}")
    # helper returns canonical repeats (no extra +1 adjustment here)
    return int(val)


def compute_id_from_aligned_read(aln) -> Tuple[int, int]:
    """Return (matches, aligned_bases) estimated from the read alignment using MD tag and CIGAR.

    We will approximate matches as aligned_bases - (mismatches + indels).
    Use get_cigar_stats from pysam if available.
    """
    # skip unmapped
    if aln.is_unmapped:
        return 0, 0
    # Use get_cigar_stats() per pysam docs. It returns two arrays; the first
    # contains counts for operations in order MIDNSHP=X and the last element is
    # the NM tag (edit distance). If no CIGAR present, arrays may be empty.
    cigstats = aln.get_cigar_stats()
    counts = cigstats[0] if cigstats is not None and len(cigstats) > 0 else []

    # counts indices: 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X,9=B,10=NM (or -1)
    # aligned on reference positions are M, =, X (indices 0,7,8)
    aligned = int(counts[0]) + int(counts[7]) + int(counts[8])

    # NM is stored as the last element; access via -1 is robust
    nm = int(counts[-1])

    # matches = aligned - nm (nm is mismatches+indels); floor at 0
    matches = max(0, aligned - nm) if aligned > 0 else 0
    return matches, aligned


def aggregate_bam(bam_path: str, tool_name: str) -> Tuple[Dict[int, List[float]], Dict[int, List[int]]]:
    """Iterate BAM and return two mappings:
    - repeats -> list of percent_id values (0-100)
    - repeats -> list of aligned lengths (int)
    """
    if pysam is None:
        raise RuntimeError('pysam is required to read BAM files')
    agg: Dict[int, List[float]] = defaultdict(list)
    len_agg: Dict[int, List[int]] = defaultdict(list)
    with pysam.AlignmentFile(bam_path, 'rb') as fh:
        for aln in fh:
            if aln.is_unmapped:
                continue
            name = aln.query_name
            try:
                if tool_name == 'tidehunter':
                    rep = get_repeats_from_tide_header(name)
                else:
                    rep = get_repeats_from_c3poa_header(name)
            except ValueError:
                # skip reads we can't parse
                continue

            matches, aligned = compute_id_from_aligned_read(aln)
            if aligned == 0:
                continue
            pct = matches / aligned * 100.0
            agg[rep].append(pct)
            len_agg[rep].append(aligned)
    return agg, len_agg


def summarize_and_write(tide_agg: Dict[int, List[float]], tide_noargs_agg: Dict[int, List[float]] | None, c3_agg: Dict[int, List[float]], tide_len: Dict[int, List[int]], tide_noargs_len: Dict[int, List[int]] | None, c3_len: Dict[int, List[int]], out_tsv: str, plot_png: str | None, length_plot: str | None = None):
    # collect all repeat keys across available datasets
    keys = set()
    keys.update(tide_agg.keys())
    keys.update(c3_agg.keys())
    if tide_noargs_agg is not None:
        keys.update(tide_noargs_agg.keys())
    keys = sorted(keys)

    with open(out_tsv, 'w') as out:
        out.write('tool\trepeats\tn_reads\tmean_id\tmedian_id\tsd_id\tmedian_aligned_len\n')
        for rep in keys:
            # tidehunter (with splint args)
            vals = tide_agg.get(rep, [])
            if vals:
                n = len(vals)
                mn = mean(vals)
                md = median(vals)
                sd = pstdev(vals) if n > 1 else 0.0
                mlen = int(median(tide_len.get(rep, [0]))) if tide_len.get(rep) else 0
                out.write(f"tidehunter\t{rep}\t{n}\t{mn:.3f}\t{md:.3f}\t{sd:.3f}\t{mlen}\n")
            # tidehunter without args (optional)
            if tide_noargs_agg is not None:
                vals = tide_noargs_agg.get(rep, [])
                if vals:
                    n = len(vals)
                    mn = mean(vals)
                    md = median(vals)
                    sd = pstdev(vals) if n > 1 else 0.0
                    mlen = int(median(tide_noargs_len.get(rep, [0]))) if tide_noargs_len and tide_noargs_len.get(rep) else 0
                    out.write(f"tidehunter_noargs\t{rep}\t{n}\t{mn:.3f}\t{md:.3f}\t{sd:.3f}\t{mlen}\n")
            # c3poa
            vals = c3_agg.get(rep, [])
            if vals:
                n = len(vals)
                mn = mean(vals)
                md = median(vals)
                sd = pstdev(vals) if n > 1 else 0.0
                mlen = int(median(c3_len.get(rep, [0]))) if c3_len.get(rep) else 0
                out.write(f"c3poa\t{rep}\t{n}\t{mn:.3f}\t{md:.3f}\t{sd:.3f}\t{mlen}\n")

    # plot median identity vs repeats if requested
    if plot_png is not None and plt is not None:
        reps = keys
        t_meds = [median(tide_agg.get(r, [0])) if tide_agg.get(r) else None for r in reps]
        c_meds = [median(c3_agg.get(r, [0])) if c3_agg.get(r) else None for r in reps]
        tn_meds = [median(tide_noargs_agg.get(r, [0])) if tide_noargs_agg and tide_noargs_agg.get(r) else None for r in reps]
        plt.figure(figsize=(8, 4))
        plt.plot(reps, [v if v is not None else float('nan') for v in t_meds], marker='o', label='tidehunter')
        if tide_noargs_agg is not None:
            plt.plot(reps, [v if v is not None else float('nan') for v in tn_meds], marker='o', label='tidehunter_noargs')
        plt.plot(reps, [v if v is not None else float('nan') for v in c_meds], marker='o', label='c3poa')
        plt.xlabel('Repeats')
        plt.ylabel('Median percent identity')
        plt.title('Median identity vs repeats')
        plt.legend()
        plt.ylim(0, 100)
        plt.tight_layout()
        plt.savefig(plot_png, dpi=150)

    # plot median aligned length vs repeats if requested (explicit path or derived from plot_png)
    if (length_plot is not None or plot_png is not None) and plt is not None:
        try:
            if length_plot is not None:
                length_plot_path = length_plot
            else:
                # derive from plot_png by appending _length.png
                if isinstance(plot_png, str) and plot_png.lower().endswith('.png'):
                    length_plot_path = plot_png[:-4] + '_length.png'
                else:
                    length_plot_path = str(plot_png) + '_length.png'

            t_len_meds = [median(tide_len.get(r, [0])) if tide_len.get(r) else None for r in reps]
            c_len_meds = [median(c3_len.get(r, [0])) if c3_len.get(r) else None for r in reps]
            tn_len_meds = [median(tide_noargs_len.get(r, [0])) if tide_noargs_len and tide_noargs_len.get(r) else None for r in reps]
            plt.figure(figsize=(8, 4))
            plt.plot(reps, [v if v is not None else float('nan') for v in t_len_meds], marker='o', label='tidehunter')
            if tide_noargs_len is not None:
                plt.plot(reps, [v if v is not None else float('nan') for v in tn_len_meds], marker='o', label='tidehunter_noargs')
            plt.plot(reps, [v if v is not None else float('nan') for v in c_len_meds], marker='o', label='c3poa')
            plt.xlabel('Repeats')
            plt.ylabel('Median aligned length (bp)')
            plt.title('Median aligned length vs repeats')
            plt.legend()
            plt.tight_layout()
            plt.savefig(length_plot_path, dpi=150)
        except Exception:
            # don't crash on plotting errors
            pass


def main(argv=None):
    args = parse_args(argv)
    tide_agg, tide_len = aggregate_bam(args.tide_bam, 'tidehunter')
    tide_noargs_agg = None
    tide_noargs_len = None
    if getattr(args, 'tide_bam_noargs', None):
        tide_noargs_agg, tide_noargs_len = aggregate_bam(args.tide_bam_noargs, 'tidehunter')
    c3_agg, c3_len = aggregate_bam(args.c3poa_bam, 'c3poa')
    summarize_and_write(tide_agg, tide_noargs_agg, c3_agg, tide_len, tide_noargs_len, c3_len, args.out_tsv, args.plot, getattr(args, 'length_plot', None))


if __name__ == '__main__':
    main()
