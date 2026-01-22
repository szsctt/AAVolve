import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from typing import Iterable, Optional

from Bio import SeqIO


def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def detect_seq_format(path: str) -> str:
    with open_maybe_gzip(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            break
    raise ValueError(f"Could not detect FASTA/FASTQ format for {path!r}")


def _clean_record(record, new_id: str):
    record.id = new_id
    record.name = new_id
    record.description = ""
    return record


def write_reference(
    *,
    reference_path: str,
    output_fasta_path: str,
) -> int:
    with open_maybe_gzip(reference_path, "rt") as ref_handle:
        ref_records = list(SeqIO.parse(ref_handle, "fasta"))
    if not ref_records:
        raise ValueError(f"No reference sequences found in {reference_path!r}")

    with open(output_fasta_path, "wt") as out_handle:
        for record in ref_records:
            record_id = record.id
            SeqIO.write(_clean_record(record, f"ref__{record_id}"), out_handle, "fasta")

    return len(ref_records)


def write_reads(
    *,
    reads_path: str,
    output_fasta_path: str,
    max_reads: int = 200,
) -> int:
    if max_reads < 0:
        raise ValueError(f"max_reads must be >= 0 (got {max_reads})")

    reads_format = detect_seq_format(reads_path)
    total_written = 0
    with open(output_fasta_path, "wt") as out_handle:
        with open_maybe_gzip(reads_path, "rt") as reads_handle:
            for index, record in enumerate(SeqIO.parse(reads_handle, reads_format)):
                if index >= max_reads:
                    break
                record_id = record.id
                SeqIO.write(
                    _clean_record(record, f"read{index+1:03d}__{record_id}"),
                    out_handle,
                    "fasta",
                )
                total_written += 1
    return total_written


def write_reference_and_reads(
    *,
    reference_path: str,
    reads_path: str,
    output_fasta_path: str,
    max_reads: int = 200,
) -> int:
    if max_reads < 0:
        raise ValueError(f"max_reads must be >= 0 (got {max_reads})")

    trimmed_format = detect_seq_format(reads_path)

    with open_maybe_gzip(reference_path, "rt") as ref_handle:
        ref_records = list(SeqIO.parse(ref_handle, "fasta"))
    if not ref_records:
        raise ValueError(f"No reference sequences found in {reference_path!r}")

    total_written = 0
    with open(output_fasta_path, "wt") as out_handle:
        for record in ref_records:
            SeqIO.write(_clean_record(record, f"ref__{record.id}"), out_handle, "fasta")
            total_written += 1

        with open_maybe_gzip(reads_path, "rt") as trimmed_handle:
            for index, record in enumerate(SeqIO.parse(trimmed_handle, trimmed_format)):
                if index >= max_reads:
                    break
                SeqIO.write(
                    _clean_record(record, f"read{index+1:03d}__{record.id}"),
                    out_handle,
                    "fasta",
                )
                total_written += 1

    return total_written


def run_mafft(*, input_fasta: str, output_fasta: str, threads: int = 1) -> None:
    if threads < 1:
        raise ValueError(f"threads must be >= 1 (got {threads})")

    command = ["mafft", "--auto", "--thread", str(threads), input_fasta]
    with open(output_fasta, "wt") as out_handle:
        subprocess.run(command, check=True, stdout=out_handle)


def run_mafft_add(
    *,
    reference_alignment_fasta: str,
    reads_fasta: str,
    output_fasta: str,
    threads: int = 1,
) -> None:
    if threads < 1:
        raise ValueError(f"threads must be >= 1 (got {threads})")

    command = ["mafft", "--thread", str(threads), "--add", reads_fasta, reference_alignment_fasta]
    with open(output_fasta, "wt") as out_handle:
        subprocess.run(command, check=True, stdout=out_handle)


def build_msa(
    *,
    reference_path: str,
    reads_path: str,
    output_msa_path: str,
    max_reads: int = 200,
    threads: int = 1,
) -> None:
    os.makedirs(os.path.dirname(output_msa_path) or ".", exist_ok=True)

    with tempfile.NamedTemporaryFile(mode="wt", suffix=".fasta", delete=False) as reference_handle:
        reference_fasta = reference_handle.name
    with tempfile.NamedTemporaryFile(mode="wt", suffix=".fasta", delete=False) as reads_handle:
        reads_fasta = reads_handle.name
    aligned_reference_fasta = None

    try:
        reference_count = write_reference(reference_path=reference_path, output_fasta_path=reference_fasta)
        reads_count = write_reads(reads_path=reads_path, output_fasta_path=reads_fasta, max_reads=max_reads)

        if reads_count == 0:
            run_mafft(input_fasta=reference_fasta, output_fasta=output_msa_path, threads=threads)
            return

        if reference_count > 1:
            with tempfile.NamedTemporaryFile(mode="wt", suffix=".fasta", delete=False) as aligned_reference_handle:
                aligned_reference_fasta = aligned_reference_handle.name
            run_mafft(input_fasta=reference_fasta, output_fasta=aligned_reference_fasta, threads=threads)

        run_mafft_add(
            reference_alignment_fasta=aligned_reference_fasta or reference_fasta,
            reads_fasta=reads_fasta,
            output_fasta=output_msa_path,
            threads=threads,
        )
    finally:
        for path in (reference_fasta, reads_fasta, aligned_reference_fasta):
            if not path:
                continue
            try:
                os.remove(path)
            except OSError:
                pass  # give up if we can't remove temp file


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a MAFFT MSA for QC: reference + first N reads from a reads file (pre- or post-trimming)."
    )
    parser.add_argument("--reference", required=True, help="Reference FASTA (optionally gzipped).")
    reads_group = parser.add_mutually_exclusive_group(required=True)
    reads_group.add_argument("--reads", help="Reads FASTA/FASTQ (optionally gzipped).")
    reads_group.add_argument(
        "--trimmed",
        help="Alias for --reads (kept for backwards compatibility).",
    )
    parser.add_argument("--output", required=True, help="Output aligned FASTA path.")
    parser.add_argument("--max-seqs", type=int, default=200, help="Number of reads to include (default: 200).")
    parser.add_argument("--threads", type=int, default=1, help="Threads for MAFFT (default: 1).")
    return parser.parse_args(argv)


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(argv)
    reads_path = args.reads or args.trimmed
    build_msa(
        reference_path=args.reference,
        reads_path=reads_path,
        output_msa_path=args.output,
        max_reads=args.max_seqs,
        threads=args.threads,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
