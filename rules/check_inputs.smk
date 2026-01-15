import os



rule validate_inputs:
    """Validate parent/reference FASTA inputs for each unique input pair."""

    input:
        parents=lambda wildcards: input_validation_map[wildcards.pair_id][0],
        reference=lambda wildcards: input_validation_map[wildcards.pair_id][1],
    output:
        ok="out/qc/input-checks/{pair_id}.txt",
    log:
        "logs/validate_inputs/{pair_id}.log"
    params:
        samples=lambda wildcards: ",".join(sorted(input_validation_samples[wildcards.pair_id])),
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        set -euo pipefail

        if ! python3 -m aavolve.check_inputs \
            --parents {input.parents} \
            --reference {input.reference} \
            --output {output.ok} \
            2> {log}; then
            echo "ERROR: input validation failed for input pair {wildcards.pair_id}" >&2
            echo "  samples: {params.samples}" >&2
            echo "  parents: {input.parents}" >&2
            echo "  reference: {input.reference}" >&2
            echo "---- validation error ----" >&2
            cat {log} >&2
            echo "--------------------------" >&2
            exit 1
        fi

        # If we got here, validation succeeded; keep a minimal success marker in the log.
        echo "OK: inputs validated for pair {wildcards.pair_id}" >> {log}
        """


rule inputs_ok:
    """Aggregate input validation and provide a single gating marker.

    Downstream rules can depend on this marker to ensure that no work starts
    until all parent/reference FASTA pairs have been validated.
    """

    input:
        input_validation_targets
    output:
        ok="out/qc/input-checks/_all.ok"
    run:
        from pathlib import Path

        out_path = Path(output.ok)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("ok\n")
