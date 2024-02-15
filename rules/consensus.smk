from scripts.snakemake_helpers import get_column_by_sample

rule consensus:
   input:
       reads = lambda wildcards: get_column_by_sample(wildcards, samples, "read_file"),
       splint = lambda wildcards: get_column_by_sample(wildcards, samples, "splint_file"),
   output:
       consensus_reads = "out/c3poa/{sample}/split/R2C2_Consensus.fasta.gz"
   params:
       dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output.consensus_reads)) + '/'
   container: "docker://szsctt/lr_c3poa"
   threads: 8
   shell:
       """
       echo "running C3POa"
       rm -r {params.dir}
       python3 /C3POa/C3POa.py \
        -r {input.reads} \
        -s {input.splint} \
        -o {params.dir} \
        -n {threads}
      
      # compress reads
      echo "compressing outputs"
      pigz -p {threads} {params.dir}/split/*

      echo "cleaning up temp files"
      # clean up temp files left behind
      rm -r {params.dir}/tmp
     """

rule filter_consensus:
    input:
        fasta = "out/c3poa/{sample}/split/R2C2_Consensus.fasta.gz"
    output:
        filt = "out/c3poa_filt/{sample}.fasta.gz"
    params:
        n_filt = lambda wildcards: int(get_column_by_sample(wildcards, samples, "min_reps"))
    container: "docker://szsctt/lr_pybio:py310"
    shell:
        """
        python3 -m scripts.filter_consensus_by_repeats -i {input.fasta} -o {output.filt} --min-repeats {params.n_filt}
        """
