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

rule consensus:
   input:
       inputs_ok = "out/qc/input-checks/_all.ok",
       reads = lambda wildcards: get_column_by_sample(wildcards, samples, "read_file"),
       splint = "out/splint/{sample}/splint.fasta",
   output:
       consensus_reads = "out/c3poa/{sample}/splint/R2C2_Consensus.fasta.gz"
   params:
       dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output.consensus_reads)) + '/'
   log:
       "logs/consensus/{sample}.log"
   container: "docker://szsctt/lr_c3poa"
   threads: 8
   shell:
       """
       echo "running C3POa"
       rm -rf {params.dir}
       python3 /C3POa/C3POa.py \
        -r {input.reads} \
        -s {input.splint} \
        -o {params.dir} \
        -n {threads}
      
      # compress reads
      echo "compressing outputs"
      pigz -p {threads} {params.dir}/splint/*

      echo "cleaning up temp files"
      # clean up temp files left behind
      rm -r {params.dir}/tmp
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

rule count_repeats:
    input:
        inputs_ok = "out/qc/input-checks/_all.ok",
        fasta = "out/c3poa/{sample}/splint/R2C2_Consensus.fasta.gz"
    output:
        counts = "out/c3poa/{sample}/repeat_counts.tsv"
    log:
        "logs/count_repeats/{sample}.log"
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m aavolve.count_RCA_repeats -i {input.fasta} -o {output.counts}
        """
