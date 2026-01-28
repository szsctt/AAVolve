import os

from aavolve.snakemake_helpers import get_column_by_sample

rule normalize_splint:
   input:
       splint = lambda wildcards: get_column_by_sample(wildcards, samples, "splint_file"),
   output:
       splint = "out/splint/{sample}/splint.fasta"
   container: "docker://szsctt/lr_pybio:py310"
   shell:
       """
       python3 -m aavolve.normalize_splint_fasta -i {input.splint} -o {output.splint}
       """

rule normalize_splint_by_input:
   input:
       inputs_ok="out/qc/input-checks/_all.ok",
       splint = lambda wildcards: npcc_cc_id_to_splint[wildcards.cc_id],
   output:
       splint = out/splint/by_input/{cc_id}/splint.fasta"
   log:
       "logs/normalize_splint_by_input/{cc_id}.log"
   container: "docker://szsctt/lr_pybio:py310"
   shell:
       """
       set -euo pipefail
       mkdir -p $(dirname {output.splint})
       python3 -m aavolve.normalize_splint_fasta -i {input.splint} -o {output.splint} 2> {log}
       """

rule consensus_by_input:
   input:
       inputs_ok = "out/qc/input-checks/_all.ok",
       reads = lambda wildcards: npcc_cc_id_to_reads[wildcards.cc_id],
       splint = "out/splint/by_input/{cc_id}/splint.fasta",
   output:
       consensus_reads = "out/c3poa/by_input/{cc_id}/splint/R2C2_Consensus.fasta.gz"
   params:
       dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output.consensus_reads)) + '/',
   log:
       "logs/consensus_by_input/{cc_id}.log"
   container: "docker://szsctt/lr_c3poa"
   threads: 8
   shell:
       """
       set -euo pipefail
       echo "running C3POa"
       rm -rf {params.dir}
       python3 /C3POa/C3POa.py \
        -r {input.reads} \
        -s {input.splint} \
        -o {params.dir} \
        -n {threads} \
        2> {log}
      
      # compress reads
      echo "compressing outputs"
      pigz -p {threads} {params.dir}/splint/*

      echo "cleaning up temp files"
      # clean up temp files left behind
      rm -r {params.dir}/tmp
      """

rule consensus:
   input:
       inputs_ok = "out/qc/input-checks/_all.ok",
       consensus_reads = lambda wildcards: f"out/c3poa/by_input/{npcc_sample_to_cc_id[wildcards.sample]}/splint/R2C2_Consensus.fasta.gz",
   output:
       consensus_reads = "out/c3poa/{sample}/splint/R2C2_Consensus.fasta.gz"
   log:
       "logs/consensus/{sample}.log"
   container: "docker://szsctt/lr_pybio:py310"
   shell:
       """
       set -euo pipefail
       mkdir -p $(dirname {output.consensus_reads})
       ln -sf $(realpath {input.consensus_reads}) {output.consensus_reads}
       """

rule filter_consensus:
    input:
        inputs_ok = "out/qc/input-checks/_all.ok",
        fasta = "out/c3poa/{sample}/splint/R2C2_Consensus.fasta.gz"
    output:
        filt = "out/c3poa_filt/{sample}.fasta.gz"
    params:
        n_filt = lambda wildcards: int(get_column_by_sample(wildcards, samples, "min_reps"))
    log:
        "logs/filter_consensus/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.filter_consensus_by_repeats -i {input.fasta} -o {output.filt} --min-repeats {params.n_filt}
        """

rule count_repeats_by_input:
    input:
        inputs_ok = "out/qc/input-checks/_all.ok",
        fasta = "out/c3poa/by_input/{cc_id}/splint/R2C2_Consensus.fasta.gz",
    output:
        counts = "out/c3poa/by_input/{cc_id}/repeat_counts.tsv",
    log:
        "logs/count_repeats_by_input/{cc_id}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        set -euo pipefail
        python3 -m aavolve.count_RCA_repeats -i {input.fasta} -o {output.counts} 2> {log}
        """

rule count_repeats:
    input:
        inputs_ok = "out/qc/input-checks/_all.ok",
        counts = lambda wildcards: f"out/c3poa/by_input/{npcc_sample_to_cc_id[wildcards.sample]}/repeat_counts.tsv",
    output:
        counts = "out/c3poa/{sample}/repeat_counts.tsv"
    log:
        "logs/count_repeats/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.counts})
        ln -sf $(realpath {input.counts}) {output.counts}
        """
