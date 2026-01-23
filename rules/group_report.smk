import os

from aavolve.snakemake_helpers import get_column_by_sample


def _samples_for_pair(pair_id: str):
    return sorted(input_validation_samples[pair_id])


def _trimmed_samples_for_pair(pair_id: str):
    trimmed = []
    for sample in _samples_for_pair(pair_id):
        wc = type("WC", (), {"sample": sample})
        if "trim" in samples.columns and bool(get_column_by_sample(wc, samples, "trim")):
            trimmed.append(sample)
    return trimmed


rule group_manifest:
    input:
        # iputs are unused, but ensure that manifest is re-run if samples change
        read_counts=lambda wildcards: [
            f"out/qc/{s}_read-counts.tsv" for s in _samples_for_pair(wildcards.pair_id)
        ],
    output:
        manifest="out/reports/group_reports/{pair_id}_manifest.tsv",
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        set -euo pipefail
        python3 -m aavolve.write_group_manifest \
            --samples-csv {input.samples_csv} \
            --pair-id {wildcards.pair_id} \
            --output {output.manifest}
        """


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
