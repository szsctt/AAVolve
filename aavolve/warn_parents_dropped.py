#!/usr/bin/env python3

import argparse
import os
import sys
from dataclasses import dataclass
from typing import Optional

import pysam

from aavolve.utils import use_open


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Warn when one or more parent sequences are absent from downstream analysis because they do not "
            "cover the reference end-to-end in the alignment."
        )
    )
    parser.add_argument("--parents-fasta", required=True, help="Parent FASTA file (multi-FASTA).")
    parser.add_argument("--bam", required=True, help="Alignment BAM of the parent FASTA against the reference (indexed).")
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA used for alignment (to get length).")
    parser.add_argument("--output", "-o", required=True, help="Output warning text file (empty if no warning).")
    return parser.parse_args(argv)


def _fasta_names(path: str) -> list[str]:
    names: list[str] = []
    with use_open(path, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            name = line[1:].strip().split()[0]
            if name:
                names.append(name)
    return names


def _single_fasta_length(path: str) -> int:
    """
    Gets length of reference from single-entry fasta file
    
    :param path: Path to fasta
    :type path: str
    :return: Description
    :rtype: int
    """
    length = 0
    with use_open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            length += len(line.strip())
    if length <= 0:
        raise ValueError(f"Reference FASTA {path!r} appears to be empty.")
    return length


def _is_primary(read: pysam.AlignedSegment) -> bool:
    """
    Return True if read is primary alignment
    
    :param read: pysam.AlignedSegment in question
    :type read: pysam.AlignedSegment
    :return: True/False value indicating if read is primary alignment
    :rtype: bool
    """
    return read.is_mapped and not read.is_secondary and not read.is_supplementary


@dataclass(frozen=True)
class _BestAlignment:
    """
    Store current best alignment
    """
    start: int
    end: int
    mapq: int
    aligned_ref_bases: int


def _best_alignment_by_qname(bam: pysam.AlignmentFile) -> dict[str, _BestAlignment]:
    """
    Get best alignment for read by qname
    
    :param bam: Description
    :type bam: pysam.AlignmentFile
    :return: Description
    :rtype: dict[str, _BestAlignment]
    """
    best: dict[str, _BestAlignment] = {}
    for read in bam.fetch(until_eof=True):
        if not _is_primary(read):
            continue
        if read.reference_start is None or read.reference_end is None:
            continue
        qname = read.query_name
        if not qname:
            continue
        # create candidate alignment
        candidate = _BestAlignment(
            start=int(read.reference_start),
            end=int(read.reference_end),
            mapq=int(read.mapping_quality or 0),
            aligned_ref_bases=int(read.reference_end - read.reference_start),
        )
        # compare previous best with candidate
        prev = best.get(qname)
        if prev is None:
            best[qname] = candidate
            continue
        if candidate.aligned_ref_bases > prev.aligned_ref_bases:
            best[qname] = candidate
        elif candidate.aligned_ref_bases == prev.aligned_ref_bases and candidate.mapq > prev.mapq:
            best[qname] = candidate
    return best


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # get names of parents expected in bam 
    parent_names = _fasta_names(args.parents_fasta)
    if not parent_names:
        raise ValueError(f"No parent sequences found in {args.parents_fasta!r}")

    # read alignment file
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        # get reference information
        if not bam.references:
            raise ValueError(f"No references found in BAM {args.bam!r}")
        ref_name = bam.references[0]
        ref_len = bam.get_reference_length(ref_name)
        if ref_len <= 0:
            ref_len = _single_fasta_length(args.reference_fasta)

        # get best alignment for each read in bam
        best_by_qname = _best_alignment_by_qname(bam)

    # check if any parents were dropped because they weren't aligned end-to-end
    dropped: list[str] = []
    for name in parent_names:
        aln = best_by_qname.get(name)
        if aln is None:
            dropped.append(name)
            continue
        if aln.start > 0 or aln.end < ref_len:
            dropped.append(name)

    with open(args.output, "wt") as handle:
        if not dropped:
            handle.write("")
            return 0
        handle.write(
            "One or more parent sequences did not align end-to-end to the reference and were dropped from "
            "parent-variant discovery.\n"
        )
        handle.write(f"- Total parents in FASTA: {len(parent_names)}\n")
        handle.write(f"- Parents aligned end-to-end: {len(parent_names) - len(dropped)}\n")
        handle.write(f"- Dropped parents ({len(dropped)}): {', '.join(dropped)}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

