import os

from aavolve.snakemake_helpers import get_column_by_sample


def _samples_for_pair(pair_id: str):
    return sorted(input_validation_samples[pair_id])


def _manifest_rows(pair_id: str):
    parent_file, reference_file = input_validation_map[pair_id]
    rows = []
    for sample in _samples_for_pair(pair_id):
        wc = type("WC", (), {"sample": sample})
        rows.append(
            {
                "sample": sample,
                "seq_tech": get_column_by_sample(wc, samples, "seq_tech"),
                "parent_file": parent_file,
                "reference_file": reference_file,
                "read_counts": f"out/qc/{sample}_read-counts.tsv",
                "assigned_parents": f"out/parents/counts/{sample}_parent-counts.tsv.gz",
                "parent_frequencies": f"out/parents/freqs/{sample}_assigned-parents_freq.tsv.gz",
                "dmat_nt_first": f"out/corrected/dmat/{sample}_first_nt-seq.tsv.gz",
                "dmat_aa_first": f"out/corrected/dmat/{sample}_first_aa-seq.tsv.gz",
                "dmat_nt_random": f"out/corrected/dmat/{sample}_random_nt-seq.tsv.gz",
                "dmat_aa_random": f"out/corrected/dmat/{sample}_random_aa-seq.tsv.gz",
            }
        )
    return rows


rule group_manifest:
    output:
        manifest="out/qc/group_reports/{pair_id}_manifest.tsv",
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
        manifest="out/qc/group_reports/{pair_id}_manifest.tsv",
        # ensure per-sample outputs exist before attempting to render
        read_counts=lambda wildcards: [
            f"out/qc/{s}_read-counts.tsv" for s in _samples_for_pair(wildcards.pair_id)
        ],
        assigned_parents=lambda wildcards: [
            f"out/parents/counts/{s}_parent-counts.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
        ],
        parent_frequencies=lambda wildcards: [
            f"out/parents/freqs/{s}_assigned-parents_freq.tsv.gz" for s in _samples_for_pair(wildcards.pair_id)
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
        report="out/qc/group_reports/{pair_id}_report.html",
        tmp_notebook="out/qc/group_reports/{pair_id}_report.ipynb",
    log:
        "logs/report_groups/{pair_id}.log"
    container: "docker://szsctt/lr_pybio:py310"
    params:
        report_basename=lambda wildcards, output: output.tmp_notebook.split('/')[-1],
    shell:
        """
        set -euo pipefail

        mkdir -p out/qc/group_reports

        papermill {input.report_template} {output.tmp_notebook} \
            -p group_id {wildcards.pair_id} \
            -p manifest {input.manifest}

        cd out/qc/group_reports
        quarto render {params.report_basename}
        """
