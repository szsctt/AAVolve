"""Input validation for AAVolve parent/reference FASTA files.

Checks performed:
  - Parents FASTA exists and parses as FASTA (supports `.gz`).
  - Reference FASTA exists and parses as FASTA (supports `.gz`).
  - Parents FASTA contains at least one record.
  - Parent sequence IDs are non-empty and unique.
  - Parent headers contain no whitespace (Biopython truncates IDs at whitespace) and no commas.
  - Parent sequences are non-empty and contain only IUPAC DNA bases: `ACGTRYSWKMBDHVN`.
  - Parent sequences are not duplicated under different IDs (identical normalized sequence).
  - Reference FASTA contains exactly one record.
  - Reference sequence is non-empty and contains only IUPAC DNA bases: `ACGTRYSWKMBDHVN`.
  - If the reference ID is also present in the parents FASTA, the sequences must match.
"""

import argparse
from pathlib import Path
import sys

from Bio import SeqIO

from aavolve.utils import use_open

IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def _normalize_seq(seq: str) -> str:
    return str(seq).upper().replace("\n", "").replace("\r", "")


def validate_inputs(parents_fasta: str, reference_fasta: str) -> None:
    """Validate parent/reference FASTA inputs.

    Raises:
        ValueError: if any input validation fails.
    """

    # paths to check
    parents_path = Path(parents_fasta)
    reference_path = Path(reference_fasta)

    # check parents and reference files exist
    if not parents_path.exists():
        raise ValueError(f"Parents FASTA not found: {parents_fasta}")
    if not reference_path.exists():
        raise ValueError(f"Reference FASTA not found: {reference_fasta}")

    # load parent records
    with use_open(parents_fasta, "rt") as handle:
        parent_records = list(SeqIO.parse(handle, "fasta"))

    # check for non-empty parents file
    if not parent_records:
        raise ValueError(f"Parents FASTA is empty: {parents_fasta}")

    # Check parents fasta names are non-empty
    parent_ids = [rec.id for rec in parent_records]
    if any(not seq_id for seq_id in parent_ids):
        raise ValueError("Parents FASTA contains an empty sequence name")

    # Check for duplicate sequence names in parent fasta
    id_counts = {}
    for seq_id in parent_ids:
        id_counts[seq_id] = id_counts.get(seq_id, 0) + 1
    dup_ids = sorted([seq_id for seq_id, count in id_counts.items() if count > 1])
    if dup_ids:
        raise ValueError(f"Duplicate sequence name(s) in parents FASTA: {', '.join(dup_ids)}")

    # Biopython's FASTA parser splits the header on whitespace and stores only
    # the first token in `record.id`. If users provide headers with whitespace,
    # that information is silently truncated. We treat that as invalid to avoid
    # confusing downstream outputs (e.g. plots/labels).
    bad_ids = [
        rec.id
        for rec in parent_records
        if rec.description != rec.id or any(c.isspace() for c in rec.id) or "," in rec.id
    ]
    if bad_ids:
        raise ValueError(
            "Parent sequence name(s) contain whitespace or commas (can break downstream plots): "
            + ", ".join(sorted(set(bad_ids)))
        )


    normalized_parent_seqs = {}
    for rec in parent_records:
        seq = _normalize_seq(rec.seq)
        if not seq:
            raise ValueError(f"Parent sequence {rec.id} is empty")
        bad_chars = sorted(set(seq) - IUPAC_DNA)
        if bad_chars:
            raise ValueError(f"Parent sequence {rec.id} has invalid base(s): {''.join(bad_chars)}")
        normalized_parent_seqs[rec.id] = seq

    # Catch the common case where the same parent sequence appears multiple times
    # under different IDs.
    seq_to_ids = {}
    for seq_id, seq in normalized_parent_seqs.items():
        seq_to_ids.setdefault(seq, []).append(seq_id)
    duplicate_sequences = [ids for ids in seq_to_ids.values() if len(ids) > 1]
    if duplicate_sequences:
        groups = ["/".join(sorted(ids)) for ids in duplicate_sequences]
        raise ValueError(
            "Duplicate parent sequences detected (different IDs, identical sequence): "
            + "; ".join(sorted(groups))
        )


    # load reference fasta
    with use_open(reference_fasta, "rt") as handle:
        reference_records = list(SeqIO.parse(handle, "fasta"))

    # check reference fasta non-empty
    if not reference_records:
        raise ValueError(f"Reference FASTA is empty: {reference_fasta}")

    # check reference contains one sequence 
    if len(reference_records) != 1:
        raise ValueError(
            f"Reference FASTA must contain exactly one sequence; found {len(reference_records)}"
        )

    # check reference sequence is non-empty
    reference_record = reference_records[0]
    reference_id = reference_record.id
    reference_seq = _normalize_seq(reference_record.seq)
    if not reference_seq:
        raise ValueError("Reference sequence is empty")

    # check for invalid bases in reference
    bad_chars = sorted(set(reference_seq) - IUPAC_DNA)
    if bad_chars:
        raise ValueError(f"Reference sequence has invalid base(s): {''.join(bad_chars)}")

    # If the reference is also one of the parents, the sequences must match.
    if reference_id in normalized_parent_seqs and normalized_parent_seqs[reference_id] != reference_seq:
        raise ValueError(
            f"Reference sequence name '{reference_id}' is present in parents FASTA, "
            "but the sequences differ"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate AAVolve parent/reference FASTA inputs")
    parser.add_argument("--parents", required=True, help="Parent sequences FASTA(.gz)")
    parser.add_argument("--reference", required=True, help="Reference sequence FASTA(.gz)")
    parser.add_argument("--output", required=True, help="Dummy output file written on success")

    args = parser.parse_args()

    try:
        validate_inputs(args.parents, args.reference)
    except ValueError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        raise SystemExit(2)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("ok\n")


if __name__ == "__main__":
    main()
