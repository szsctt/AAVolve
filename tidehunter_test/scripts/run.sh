#!/bin/bash
set -euo pipefail


# Check if project is run from project root or not - if not, change to project root
if [ ! -f "tidehunter_test/scripts/run.sh" ]; then
  echo "Please run this script from the project root directory."
  exit 1
fi

DATADIR="tidehunter_test/data"
READS=$DATADIR/SRR29543177.fastq
AAV2_REF="tidehunter_test/data/aav2.fa"

# create micromamba env
eval "$(micromamba shell hook -s bash)"
micromamba create -n tidehunter_env -c conda-forge -c bioconda tidehunter wget sra-tools minimap2 samtools pysam matplotlib -y

# Create data directory if it doesn't exist
mkdir -p $DATADIR

# Download data from SRA - AAV2 sample with R2C2 sequencing
if [ ! -f "$READS" ]; then
    echo "Downloading data from SRA..."
    micromamba run -n tidehunter_env fastq-dump SRR29543177
    mv SRR29543177.fastq $READS
fi

# Run TideHunter on the downloaded data
SPLINT=$DATADIR/splint.fa
OUTTIDE=tidehunter_test/output/tidehunter
mkdir -p $OUTTIDE
/usr/bin/time -v micromamba run -n tidehunter_env \
    TideHunter \
    -5 $SPLINT \
    -3 $SPLINT \
    -t 20 \
    --longest \
    -F \
    $READS  > $OUTTIDE/tidehunter_output.fa 2> $OUTTIDE/tidehunter_time.txt

# Run TideHunter without -5/-3 (no splint arguments)
# output FASTA and timing/logs
/usr/bin/time -v micromamba run -n tidehunter_env \
    TideHunter \
    -t 20 \
    --longest \
    $READS > $OUTTIDE/tidehunter_output_noargs.fa 2> $OUTTIDE/tidehunter_noargs_time.txt


# Run c3poa 
OUTC3POA=tidehunter_test/output/c3poa/
mkdir -p $OUTC3POA
if [ ! -f lr_c3poa_latest.sif ]; then
    echo "Downloading C3POa Singularity image..."
    singularity pull docker://szsctt/lr_c3poa:latest
fi
#/usr/bin/time -v singularity exec lr_c3poa_latest.sif \
#    /C3POa/C3POa.py \
#    -r $READS \
#    -s $SPLINT \
#    -o $OUTC3POA \
#    -n 20 > $OUTC3POA/c3poa_time.txt 2>&1


# plot repeats
mkdir -p tidehunter_test/output/plots
python3 tidehunter_test/scripts/plot_repeats.py \
  --tidehunter tidehunter_test/output/tidehunter/tidehunter_output.fa \
    --tidehunter-noargs tidehunter_test/output/tidehunter/tidehunter_output_noargs.fa \
  --c3poa tidehunter_test/output/c3poa/split/R2C2_Consensus.fasta \
  --out tidehunter_test/output/plots/repeats_compare.png \
  --log-x

# align both to aav2 reference
ALNTIDE=tidehunter_test/output/alignments/tidehunter_vs_aav2.bam
ALNC3POA=tidehunter_test/output/alignments/c3poa_vs_aav2.bam
ALNTIDE_NOARGS=tidehunter_test/output/alignments/tidehunter_noargs_vs_aav2.bam
mkdir -p tidehunter_test/output/alignments
micromamba run -n tidehunter_env \
    minimap2 -ax map-ont $AAV2_REF $OUTTIDE/tidehunter_output.fa --MD | samtools sort -o $ALNTIDE
micromamba run -n tidehunter_env samtools index $ALNTIDE
micromamba run -n tidehunter_env \
    minimap2 -ax map-ont $AAV2_REF $OUTC3POA/split/R2C2_Consensus.fasta --MD | samtools sort -o $ALNC3POA
micromamba run -n tidehunter_env samtools index $ALNC3POA
micromamba run -n tidehunter_env \
    minimap2 -ax map-ont $AAV2_REF $OUTTIDE/tidehunter_output_noargs.fa --MD | samtools sort -o $ALNTIDE_NOARGS
micromamba run -n tidehunter_env samtools index $ALNTIDE_NOARGS

micromamba run -n tidehunter_env samtools flagstats $ALNTIDE > tidehunter_test/output/alignments/tidehunter_flagstats.txt
micromamba run -n tidehunter_env samtools flagstats $ALNC3POA > tidehunter_test/output/alignments/c3poa_flagstats.txt
micromamba run -n tidehunter_env samtools flagstats $ALNTIDE_NOARGS > tidehunter_test/output/alignments/tidehunter_noargs_flagstats.txt


# compute accuracy by repeat count and write TSV + PNG
micromamba run -n tidehunter_env \
python3 tidehunter_test/scripts/compute_accuracy_by_repeats.py \
  --tide-bam $ALNTIDE \
  --tide-bam-noargs $ALNTIDE_NOARGS \
  --c3poa-bam $ALNC3POA \
  --out-tsv tidehunter_test/output/plots/accuracy_by_repeats.tsv \
  --plot tidehunter_test/output/plots/accuracy_by_repeats.png \
  --length-plot tidehunter_test/output/plots/length_by_repeats.png 




