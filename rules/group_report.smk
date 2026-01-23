import os

from aavolve.snakemake_helpers import get_column_by_sample


def _samples_for_pair(pair_id: str):
    return sorted(input_validation_samples[pair_id])


def _manifest_rows(pair_id: str):
    parent_file, reference_file = input_validation_map[pair_id]
    rows = []
    for sample in _samples_for_pair(pair_id):
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


def _trimmed_samples_for_pair(pair_id: str):
    trimmed = []
    for sample in _samples_for_pair(pair_id):
        wc = type("WC", (), {"sample": sample})
        if "trim" in samples.columns and bool(get_column_by_sample(wc, samples, "trim")):
            trimmed.append(sample)
    return trimmed


rule group_manifest:
    output:
        manifest="out/reports/group_reports/{pair_id}_manifest.tsv",
    run:
        from pathlib import Path
        rows = _manifest_rows(wildcards.pair_id)
        import csv

        Path(output.manifest).parent.mkdir(parents=True, exist_ok=True)
        with open(output.manifest, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


rule group_report:
    input:
        manifest="out/reports/group_reports/{pair_id}_manifest.tsv",
        # ensure per-sample outputs exist before attempting to render
        read_counts=lambda wildcards: [
            f"out/qc/{s}_read-counts.tsv" for s in _samples_for_pair(wildcards.pair_id)
        ],
        coverage_depth=lambda wildcards: [
            f"out/qc/coverage/{s}_depth.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        assigned_parents=lambda wildcards: [
            f"out/parents/counts/{s}_parent-counts.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        parent_frequencies=lambda wildcards: [
            f"out/parents/freqs/{s}_assigned-parents_freq.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        variant_freq_all=lambda wildcards: [
            f"out/variants/frequency/{s}_all.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        parents_dropped_warning=lambda wildcards: [
            f"out/qc/warnings/{s}_parents_dropped.txt" for s in _samples_for_pair(wildcards.pair_id)
        ],
        variant_window_warning=lambda wildcards: [
            f"out/qc/warnings/{s}_non_parental_outside_window.txt" for s in _samples_for_pair(wildcards.pair_id)
        ],
        pretrim_msa=lambda wildcards: [
            f"out/qc/trimming/{s}.pretrim.mafft.fasta" for s in _samples_for_pair(wildcards.pair_id)
        ],
        posttrim_msa=lambda wildcards: [
            f"out/qc/trimming/{s}.mafft.fasta" for s in _trimmed_samples_for_pair(wildcards.pair_id)
        ],
        dmat_nt_first=lambda wildcards: [
            f"out/corrected/dmat/{s}_first_nt-seq.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        dmat_aa_first=lambda wildcards: [
            f"out/corrected/dmat/{s}_first_aa-seq.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        dmat_nt_random=lambda wildcards: [
            f"out/corrected/dmat/{s}_random_nt-seq.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        dmat_aa_random=lambda wildcards: [
            f"out/corrected/dmat/{s}_random_aa-seq.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        report_template=lambda wildcards: os.path.join(workflow.basedir, "aavolve/group_report.ipynb"),
    output:
        report="out/reports/group_reports/{pair_id}_report.html",
        tmp_notebook=temp("out/reports/group_reports/{pair_id}_report.ipynb"),
    log:
        "logs/report_groups/{pair_id}.log"
    container: "docker://szsctt/lr_pybio:py310"
    params:
        report_basename=lambda wildcards, output: output.tmp_notebook.split('/')[-1],
        report_dir = lambda wildcards, output: os.path.dirname(output.report),
    shell:
        """
        set -euo pipefail

        mkdir -p out/qc/group_reports

        papermill {input.report_template} {output.tmp_notebook} \
            -p group_id {wildcards.pair_id} \
            -p manifest {input.manifest}

        cd {params.report_dir}
        unset QUARTO_DENO DENO
        quarto render {params.report_basename}
        """
