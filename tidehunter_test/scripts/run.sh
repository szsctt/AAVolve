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
micromamba create -n tidehunter_env -c conda-forge -c bioconda tidehunter wget sra-tools minimap2 samtools pysam matplotlib mafft -y

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
#/usr/bin/time -v micromamba run -n tidehunter_env \
    TideHunter \
    -5 $SPLINT \
    -3 $SPLINT \
    -t 20 \
    --longest \
    -F \
    $READS  > $OUTTIDE/tidehunter_output.fa 2> $OUTTIDE/tidehunter_time.txt

# Run TideHUnter with 5' and 3' sequences determined from reads
/usr/bin/time -v micromamba run -n tidehunter_env \
    TideHunter \
    -5 tidehunter_test/data/splint_5prime.fa \
    -3 tidehunter_test/data/splint_3prime.fa \
    -t 20 \
    --longest \
    -F \
    $READS  > $OUTTIDE/tidehunter_empiricalseqs_output.fa 2> $OUTTIDE/tidehunter_empiricalseqs_time.txt

# Run TideHunter without -5/-3 (no splint arguments)
# output FASTA and timing/logs
/usr/bin/time -v micromamba run -n tidehunter_env \
    TideHunter \
    -t 20 \
    --longest \
    $READS > $OUTTIDE/tidehunter_noprimers_output.fa 2> $OUTTIDE/tidehunter_noprimers_time.txt


# Run c3poa 
OUTC3POA=tidehunter_test/output/c3poa/
mkdir -p $OUTC3POA
if [ ! -f lr_c3poa_latest.sif ]; then
    echo "Downloading C3POa Singularity image..."
    singularity pull docker://szsctt/lr_c3poa:latest
fi
/usr/bin/time -v singularity exec lr_c3poa_latest.sif \
    /C3POa/C3POa.py \
    -r $READS \
    -s $SPLINT \
    -o $OUTC3POA \
    -n 20 > $OUTC3POA/c3poa_time.txt 2>&1


# plot repeats
mkdir -p tidehunter_test/output/plots
micromamba run -n tidehunter_env \
python3 tidehunter_test/scripts/plot_repeats.py \
  --tidehunter tidehunter_test/output/tidehunter/tidehunter_empiricalseqs_output.fa \
  --tidehunter-noargs tidehunter_test/output/tidehunter/tidehunter_noprimers_output.fa \
  --c3poa tidehunter_test/output/c3poa/split/R2C2_Consensus.fasta \
  --out tidehunter_test/output/plots/repeats_compare.png \
  --log-x

# align to aav2 reference (use helper to avoid repetition)
ALNTIDE=tidehunter_test/output/alignments/tidehunter_empiricalseqs_vs_aav2.bam
ALNC3POA=tidehunter_test/output/alignments/c3poa_vs_aav2.bam
ALNTIDE_noprimers=tidehunter_test/output/alignments/tidehunter_noprimers_vs_aav2.bam
mkdir -p tidehunter_test/output/alignments

# Helper: align FASTA to reference, sort to BAM, index and write flagstats
align_and_index() {
  local in_fa="$1"
  local out_bam="$2"
  local label="$3"

  if [ ! -f "$in_fa" ]; then
    echo "Warning: input FASTA not found: $in_fa -- skipping $label"
    return 0
  fi

  echo "Aligning $in_fa -> $out_bam"
  micromamba run -n tidehunter_env minimap2 -ax map-ont "$AAV2_REF" "$in_fa" --MD | samtools sort -o "$out_bam"
  micromamba run -n tidehunter_env samtools index "$out_bam"
  micromamba run -n tidehunter_env samtools flagstats "$out_bam" > "tidehunter_test/output/alignments/${label}_flagstats.txt"
}

# run alignments for each condition
align_and_index "$OUTTIDE/tidehunter_empiricalseqs_output.fa" "$ALNTIDE" tidehunter_empiricalseqs
align_and_index "$OUTC3POA/split/R2C2_Consensus.fasta" "$ALNC3POA" c3poa
align_and_index "$OUTTIDE/tidehunter_noprimers_output.fa" "$ALNTIDE_noprimers" tidehunter_noprimers


# compute accuracy by repeat count and write TSV + PNG
micromamba run -n tidehunter_env \
python3 tidehunter_test/scripts/compute_accuracy_by_repeats.py \
  --tide-bam $ALNTIDE \
  --tide-bam-noargs $ALNTIDE_noprimers \
  --c3poa-bam $ALNC3POA \
  --out-tsv tidehunter_test/output/plots/accuracy_by_repeats.tsv \
  --plot tidehunter_test/output/plots/accuracy_by_repeats.png \
  --length-plot tidehunter_test/output/plots/length_by_repeats.png \
  --cover-plot tidehunter_test/output/plots/cover_by_repeats.png


# --- Create small sample FASTAs (reference + first 100 reads) and run MAFFT MSAs
mkdir -p tidehunter_test/output/msa
REF_FASTA=$DATADIR/aav2.fa
TIDE_EMP_FASTA=$OUTTIDE/tidehunter_empiricalseqs_output.fa
TIDE_FASTA=$OUTTIDE/tidehunter_output.fa
TIDE_noprimers_FASTA=$OUTTIDE/tidehunter_noprimers_output.fa
C3_FASTA=$OUTC3POA/split/R2C2_Consensus.fasta

echo "Generating sample FASTAs (reference + up to 100 reads) in tidehunter_test/output/msa"

write_sample() {
  local src_fasta="$1"; local out_fasta="$2"
  if [ ! -f "$src_fasta" ]; then
    echo "  source FASTA missing: $src_fasta -- skipping"
    return
  fi
  # write reference first
  cp "$REF_FASTA" "$out_fasta"
  # append first 100 records from source (awk splits on '>')
  awk 'BEGIN{RS=">"; ORS=""} NR>1 && NR<=101{print ">"$0}' "$src_fasta" >> "$out_fasta"
  echo "  wrote $out_fasta"
}

write_sample "$TIDE_FASTA" tidehunter_test/output/msa/tidehunter_sample.fa
write_sample "$TIDE_EMP_FASTA" tidehunter_test/output/msa/tidehunter_empiricalseqs_sample.fa
write_sample "$TIDE_noprimers_FASTA" tidehunter_test/output/msa/tidehunter_noprimers_sample.fa
write_sample "$C3_FASTA" tidehunter_test/output/msa/c3poa_sample.fa

echo "Attempting to run MAFFT on samples (will skip if mafft not found)"
run_mafft() {
  local in_fa="$1"; local out_aln="$2"
    micromamba run -n tidehunter_env mafft --auto "$in_fa" > "$out_aln"
}

run_mafft tidehunter_test/output/msa/tidehunter_sample.fa tidehunter_test/output/msa/tidehunter_sample.aln.fa
run_mafft tidehunter_test/output/msa/tidehunter_empiricalseqs_sample.fa tidehunter_test/output/msa/tidehunter_empiricalseqs_sample.aln.fa
run_mafft tidehunter_test/output/msa/tidehunter_noprimers_sample.fa tidehunter_test/output/msa/tidehunter_noprimers_sample.aln.fa
run_mafft tidehunter_test/output/msa/c3poa_sample.fa tidehunter_test/output/msa/c3poa_sample.aln.fa





