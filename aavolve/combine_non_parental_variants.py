import argparse
import os
import sys
from dataclasses import dataclass
from typing import Iterable, Optional

from aavolve.utils import get_header, get_reference_name, get_variant, read_variant_file, use_open


@dataclass(frozen=True)
class VariantRecord:
    variant: object
    aa_change: bool


def _sort_key(variant) -> tuple[int, int, str, str]:
    var_id = variant.var_id()
    pos_text, vartype = var_id.split(":")
    if "_" in pos_text:
        start_text, end_text = pos_text.split("_", 1)
        start = int(start_text)
        end = int(end_text)
    else:
        start = int(pos_text)
        end = start
    return (start, end, vartype, str(variant))


def combine_non_parental_variants(
    input_files: Iterable[str],
    output_file: str,
) -> int:
    input_files = [str(path) for path in input_files]
    if not input_files:
        raise ValueError("Expected at least one input file")

    header = get_header(input_files[0]).rstrip("\n")
    columns = header.split("\t")
    long_format = columns and columns[0] == "reference_name"

    reference_name = None
    if long_format:
        reference_name = get_reference_name(input_files[0], shorter_behaviour="error")

    variants: dict[str, VariantRecord] = {}
    for input_file in input_files:
        for row in read_variant_file(input_file):
            variant = get_variant(row)
            variants.setdefault(
                str(variant),
                VariantRecord(variant=variant, aa_change=variant.changes_aa),
            )
            if reference_name is None and long_format:
                reference_name = row.get("reference_name")

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with use_open(output_file, "wt", newline="") as handle:
        handle.write(header + "\n")
        for index, record in enumerate(sorted(variants.values(), key=lambda r: _sort_key(r.variant))):
            name = f"non_parental_{index}"
            if long_format:
                if reference_name is None:
                    raise ValueError("reference_name is required for long-format variant outputs")
                handle.write(record.variant.print_line(query_name=name, ref_name=reference_name))
            else:
                handle.write(record.variant.print_line(query_name=name))

    return len(variants)


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine non-parental variant TSVs across samples.")
    parser.add_argument("inputs", nargs="+", help="Input TSV/TSV.GZ files of non-parental variants.")
    parser.add_argument("--output", "-o", required=True, help="Output TSV/TSV.GZ path.")
    return parser.parse_args(argv)


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(argv)
    combine_non_parental_variants(args.inputs, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

