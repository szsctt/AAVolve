from __future__ import annotations

import argparse
import csv
from pathlib import Path

from aavolve.get_samples import get_samples
from aavolve.snakemake_helpers import build_input_validation_targets, get_column_by_sample


def _manifest_rows(samples_csv: str, pair_id: str) -> list[dict[str, object]]:
    samples = get_samples({"samples": samples_csv})
    input_validation_map, input_validation_samples, _ = build_input_validation_targets(samples)
    if pair_id not in input_validation_map:
        raise ValueError(f"Unknown pair_id: {pair_id}")
    parent_file, reference_file = input_validation_map[pair_id]

    rows: list[dict[str, object]] = []
    for sample in sorted(input_validation_samples[pair_id]):
        wc = type("WC", (), {"sample": sample})
        trim_enabled = bool(get_column_by_sample(wc, samples, "trim")) if "trim" in samples.columns else False

        rows.append(
            {
                "sample": sample,
                "seq_tech": get_column_by_sample(wc, samples, "seq_tech"),
                "parent_file": parent_file,
                "reference_file": reference_file,
                "trim": get_column_by_sample(wc, samples, "trim") if "trim" in samples.columns else False,
                "adapter_5": get_column_by_sample(wc, samples, "adapter_5") if "adapter_5" in samples.columns else "",
                "adapter_3": get_column_by_sample(wc, samples, "adapter_3") if "adapter_3" in samples.columns else "",
                "anchors": get_column_by_sample(wc, samples, "anchors") if "anchors" in samples.columns else 0,
                "require_end_to_end_alignment": get_column_by_sample(wc, samples, "require_end_to_end_alignment")
                if "require_end_to_end_alignment" in samples.columns
                else False,
                "include_non_parental": get_column_by_sample(wc, samples, "include_non_parental"),
                "non_parental_freq": get_column_by_sample(wc, samples, "non_parental_freq")
                if "non_parental_freq" in samples.columns
                else "",
                "group_vars": get_column_by_sample(wc, samples, "group_vars") if "group_vars" in samples.columns else "",
                "group_vars_dist": get_column_by_sample(wc, samples, "group_vars_dist")
                if "group_vars_dist" in samples.columns
                else "",
                "max_group_distance": get_column_by_sample(wc, samples, "max_group_distance")
                if "max_group_distance" in samples.columns
                else "",
                "minimap2_params": get_column_by_sample(wc, samples, "minimap2_params")
                if "minimap2_params" in samples.columns
                else "",
                "read_counts": f"out/qc/{sample}_read-counts.tsv",
                "coverage_depth": f"out/qc/coverage/{sample}_depth.tsv.gz",
                "assigned_parents": f"out/parents/counts/{sample}_parent-counts.tsv.gz",
                "parent_frequencies": f"out/parents/freqs/{sample}_assigned-parents_freq.tsv.gz",
                "variant_freq_all": f"out/variants/frequency/{sample}_all.tsv.gz",
                "parents_dropped_warning": f"out/qc/warnings/{sample}_parents_dropped.txt",
                "variant_window_warning": f"out/qc/warnings/{sample}_non_parental_outside_window.txt",
                "pretrim_msa": f"out/qc/trimming/{sample}.pretrim.mafft.fasta",
                "posttrim_msa": f"out/qc/trimming/{sample}.mafft.fasta" if trim_enabled else "",
                "dmat_nt_first": f"out/corrected/dmat/{sample}_first_nt-seq.tsv.gz",
                "dmat_aa_first": f"out/corrected/dmat/{sample}_first_aa-seq.tsv.gz",
                "dmat_nt_random": f"out/corrected/dmat/{sample}_random_nt-seq.tsv.gz",
                "dmat_aa_random": f"out/corrected/dmat/{sample}_random_aa-seq.tsv.gz",
            }
        )

    return rows


class _ListWriter:
    def __init__(self) -> None:
        self.parts: list[str] = []

    def write(self, s: str) -> int:
        self.parts.append(s)
        return len(s)

    def getvalue(self) -> str:
        return "".join(self.parts)


def write_manifest(samples_csv: str, pair_id: str, output_path: str) -> None:
    rows = _manifest_rows(samples_csv=samples_csv, pair_id=pair_id)

    writer_buf = _ListWriter()
    writer = csv.DictWriter(writer_buf, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    new_contents = writer_buf.getvalue()

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        old_contents = out_path.read_text()
        if old_contents == new_contents:
            return
    out_path.write_text(new_contents)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-csv", required=True)
    parser.add_argument("--pair-id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    write_manifest(samples_csv=args.samples_csv, pair_id=args.pair_id, output_path=args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
