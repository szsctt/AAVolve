#!/bin/bash
set -euo pipefail


# Check if project is run from project root or not - if not, change to project root
if [ ! -f "trimming_test/scripts/run.sh" ]; then
  echo "Please run this script from the project root directory."
  exit 1
fi

DATADIR="trimming_test/data"
READS=$DATADIR/SRR29543177.fastq
AAV2_REF="trimming_test/data/aav2.fa"

# create micromamba env
eval "$(micromamba shell hook -s bash)"
micromamba create -n trimming_env -c conda-forge -c bioconda sra-tools mafft cutadapt trimmomatic bbmap -y


# Create data directory if it doesn't exist
mkdir -p $DATADIR

# Download data from SRA - AAV2 sample with R2C2 sequencing
if [ ! -f "$READS" ]; then
    echo "Downloading data from SRA..."
    micromamba run -n trimming_env fastq-dump SRR29543177
    mv SRR29543177.fastq $READS
fi

# Run C3POa to generate consensus reads
OUTC3POA=trimming_test/output/c3poa/
SPLINT=$DATADIR/splint.fa
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
    -n 40 > $OUTC3POA/c3poa_time.txt 2>&1

READS="$OUTC3POA/splint/R2C2_Consensus.fasta"

# Convert FASTA to FASTQ (Trimmomatic requires FASTQ format)
# Add dummy quality scores (all high quality 'I' = Q40)
READS_FQ="${READS%.fasta}.fastq"
if [ ! -f "$READS_FQ" ]; then
    echo "Converting FASTA to FASTQ with dummy quality scores..."
    awk '{
        if (substr($0,1,1)==">") {
            if (seq!="") {
                print "@" name
                print seq
                print "+"
                for (i=1; i<=length(seq); i++) printf "I"
                print ""
            }
            name=substr($0,2)
            seq=""
        } else {
            seq=seq $0
        }
    }
    END {
        if (seq!="") {
            print "@" name
            print seq
            print "+"
            for (i=1; i<=length(seq); i++) printf "I"
            print ""
        }
    }' "$READS" > "$READS_FQ"
fi

# run trimmomatic with increased heap size for long reads
# The trimmomatic wrapper script accepts -Xmx arguments
OUTTRIM=trimming_test/output/trimmomatic/
ADAPTERS=$DATADIR/adapters.fa
mkdir -p $OUTTRIM
/usr/bin/time -v micromamba run -n trimming_env \
    trimmomatic -Xmx8g SE \
    -threads 4 -phred33 \
    $READS_FQ \
    $OUTTRIM/consensus_reads_trimmed.fastq \
    ILLUMINACLIP:$ADAPTERS:2:30:10 > $OUTTRIM/trimmomatic_time.txt 2>&1

# run cutadapt (works with both FASTA and FASTQ, but use FASTQ for consistency)
OUTCUTADAPT=trimming_test/output/cutadapt/
mkdir -p $OUTCUTADAPT
/usr/bin/time -v micromamba run -n trimming_env \
    cutadapt -a TTGCTTGTTAAATCA -g TGATTTAAATCAGGT \
    --revcomp --trimmed-only \
    -o $OUTCUTADAPT/consensus_reads_cutadapt.fastq \
    $READS_FQ > $OUTCUTADAPT/cutadapt_time.txt 2>&1

# run cutadapt with non-internal adapter sequences (X prefix means non-internal)
OUTCUTADAPT_NONINT=trimming_test/output/cutadapt_noninternal/
mkdir -p $OUTCUTADAPT_NONINT
/usr/bin/time -v micromamba run -n trimming_env \
    cutadapt -a TTGCTTGTTAAATCAX -g XTGATTTAAATCAGGT \
    --revcomp --trimmed-only \
    -o $OUTCUTADAPT_NONINT/consensus_reads_cutadapt_noninternal.fastq \
    $READS_FQ > $OUTCUTADAPT_NONINT/cutadapt_noninternal_time.txt 2>&1

# run cutadapt with linked adapters
OUTCUTADAPT_LINKED=trimming_test/output/cutadapt_linked/
mkdir -p $OUTCUTADAPT_LINKED
/usr/bin/time -v micromamba run -n trimming_env \
    cutadapt -a TGATTTAAATCAGGT...TTGCTTGTTAAATCA  \
    --revcomp --trimmed-only \
    -o $OUTCUTADAPT_LINKED/consensus_reads_cutadapt_linked.fastq \
    $READS_FQ > $OUTCUTADAPT_LINKED/cutadapt_linked_time.txt 2>&1

# Run bbduk to trim adapters (left trim)
OUTBBMAP=trimming_test/output/bbmap/
mkdir -p $OUTBBMAP
/usr/bin/time -v micromamba run -n trimming_env \
    bbduk.sh in=$READS_FQ out=$OUTBBMAP/consensus_reads_bbduk_ktriml.fastq \
    ref=$ADAPTERS ktrim=l k=13 mink=5 hdist=1 rcomp=t \
    tpe tbo threads=4 > $OUTBBMAP/bbduk_time_triml.txt 2>&1

# Run bbduk to trim adapters (right trim)
/usr/bin/time -v micromamba run -n trimming_env \
    bbduk.sh in=$OUTBBMAP/consensus_reads_bbduk_ktriml.fastq \
    out=$OUTBBMAP/consensus_reads_bbduk_ktrimr.fastq \
    ref=$ADAPTERS ktrim=r k=13 mink=5 hdist=1 rcomp=t \
    tpe tbo threads=4 > $OUTBBMAP/bbduk_time_trimr.txt 2>&1

# Create multiple sequence alignments for each trimming result
# Extract first 100 reads, concatenate with AAV2 reference, and align
echo "Creating multiple sequence alignments..."
MSADIR=trimming_test/output/msa
mkdir -p $MSADIR

# Function to convert FASTQ to FASTA and get first 100 reads
fastq_to_fasta_first100() {
    local input=$1
    local output=$2
    # FASTQ records are 4 lines each, so we need 400 lines for 100 reads
    # Use head first to avoid SIGPIPE from awk when processing large files
    head -n 400 "$input" | awk 'NR % 4 == 1 {gsub(/^@/, ">"); print} NR % 4 == 2 {print}' > "$output"
}

# Process each trimming tool's output
for tool_output in \
    "$OUTTRIM/consensus_reads_trimmed.fastq:trimmomatic" \
    "$OUTCUTADAPT/consensus_reads_cutadapt.fastq:cutadapt" \
    "$OUTCUTADAPT_NONINT/consensus_reads_cutadapt_noninternal.fastq:cutadapt_noninternal" \
    "$OUTCUTADAPT_LINKED/consensus_reads_cutadapt_linked.fastq:cutadapt_linked" \
    "$OUTBBMAP/consensus_reads_bbduk_ktrimr.fastq:bbduk_trimr"
do
    # Split the path and tool name
    IFS=':' read -r fastq_file tool_name <<< "$tool_output"
    
    # Skip if file doesn't exist
    if [ ! -f "$fastq_file" ]; then
        echo "Skipping $tool_name - file not found: $fastq_file"
        continue
    fi
    
    echo "Processing $tool_name..."
    
    # Extract first 100 reads and convert to FASTA
    temp_reads="$MSADIR/${tool_name}_first100.fa"
    fastq_to_fasta_first100 "$fastq_file" "$temp_reads"
    
    # Concatenate AAV2 reference with the 100 reads
    combined_fa="$MSADIR/${tool_name}_with_aav2.fa"
    cat "$AAV2_REF" "$temp_reads" > "$combined_fa"
    
    # Run MAFFT alignment
    echo "Running MAFFT alignment for $tool_name..."
    /usr/bin/time -v micromamba run -n trimming_env \
        mafft --auto  --adjustdirection --thread 4 "$combined_fa" > "$MSADIR/${tool_name}_aligned.fa" 2> "$MSADIR/${tool_name}_mafft.log"
    
    # Clean up temporary files
    rm "$temp_reads"
    
    echo "Completed alignment for $tool_name"
done

echo "All alignments completed!"
echo "Output files are in: $MSADIR"

# Summary of trimming results
echo ""
echo "========================================"
echo "TRIMMING SUMMARY"
echo "========================================"
echo ""

# Function to parse trimmomatic time file
parse_trimmomatic() {
    local file=$1
    local tool_name=$2
    
    if [ ! -f "$file" ]; then
        echo "$tool_name: File not found"
        return
    fi
    
    # Extract runtime (elapsed wall clock time)
    runtime=$(grep "Elapsed (wall clock) time" "$file" | sed 's/.*: //')
    
    # Extract read counts from Trimmomatic output
    input_reads=$(grep "Input Reads:" "$file" | sed 's/Input Reads: //' | awk '{print $1}')
    surviving_reads=$(grep "Input Reads:" "$file" | sed 's/.*Surviving: //' | awk '{print $1}')
    dropped_reads=$(grep "Input Reads:" "$file" | sed 's/.*Dropped: //' | awk '{print $1}')
    
    echo "$tool_name:"
    echo "  Runtime: $runtime"
    echo "  Total reads: $input_reads"
    echo "  Surviving reads: $surviving_reads"
    echo "  Dropped reads: $dropped_reads"
    echo ""
}

# Function to parse cutadapt time file
parse_cutadapt() {
    local file=$1
    local tool_name=$2
    
    if [ ! -f "$file" ]; then
        echo "$tool_name: File not found"
        return
    fi
    
    # Extract runtime (elapsed wall clock time)
    runtime=$(grep "Elapsed (wall clock) time" "$file" | sed 's/.*: //')
    
    # Extract read counts from cutadapt summary
    total_reads=$(grep "Total reads processed:" "$file" | sed 's/.*: *//' | sed 's/,//g' | awk '{print $1}')
    reads_with_adapters=$(grep "Reads with adapters:" "$file" | sed 's/.*: *//' | sed 's/,//g' | awk '{print $1}')
    reads_written=$(grep "Reads written (passing filters):" "$file" | sed 's/.*: *//' | sed 's/,//g' | awk '{print $1}')
    reads_trimmed=$reads_with_adapters
    
    echo "$tool_name:"
    echo "  Runtime: $runtime"
    echo "  Total reads: $total_reads"
    echo "  Surviving reads: $reads_written"
    echo "  Reads with adapters trimmed: $reads_trimmed"
    echo ""
}

# Function to parse bbduk time file
parse_bbduk() {
    local file=$1
    local tool_name=$2
    
    if [ ! -f "$file" ]; then
        echo "$tool_name: File not found"
        return
    fi
    
    # Extract runtime (elapsed wall clock time)
    runtime=$(grep "Elapsed (wall clock) time" "$file" | sed 's/.*: //')
    
    # Extract read counts from bbduk output
    input_reads=$(grep "^Input:" "$file" | awk '{print $2}')
    ktrimmed_reads=$(grep "^KTrimmed:" "$file" | awk '{print $2}')
    result_reads=$(grep "^Result:" "$file" | awk '{print $2}')
    
    echo "$tool_name:"
    echo "  Runtime: $runtime"
    echo "  Total reads: $input_reads"
    echo "  Surviving reads: $result_reads"
    echo "  Reads trimmed: $ktrimmed_reads"
    echo ""
}

# Parse all trimming tool results
parse_trimmomatic "$OUTTRIM/trimmomatic_time.txt" "Trimmomatic"
parse_cutadapt "$OUTCUTADAPT/cutadapt_time.txt" "Cutadapt (internal)"
parse_cutadapt "$OUTCUTADAPT_NONINT/cutadapt_noninternal_time.txt" "Cutadapt (non-internal)"
parse_cutadapt "$OUTCUTADAPT_LINKED/cutadapt_linked_time.txt" "Cutadapt (linked)"
parse_bbduk "$OUTBBMAP/bbduk_time_triml.txt" "BBDuk (left trim)"
parse_bbduk "$OUTBBMAP/bbduk_time_trimr.txt" "BBDuk (right trim)"

echo "========================================"