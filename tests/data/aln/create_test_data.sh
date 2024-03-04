 AAV2='../references/wtAAV2.fa'

  minimap2 -ax map-hifi -B 1.5 --MD $AAV2 $AAV2| samtools sort -o aav2.bam -
 samtools index aav2.bam
 
 minimap2 -ax map-hifi -B 1.5 --MD --end-bonus 5 $AAV2 ../references/wtAAV2_N496D.fa | samtools sort -o aav2_N496D.bam -
 samtools index aav2.bam

 minimap2 -ax map-hifi -B 1.5 --MD --end-bonus 5 $AAV2 pb-shuf.fq | samtools sort -o pb-shuf.bam -
 samtools index pb-shuf.bam

# manually edit the sam file to add secondary and supplementary alignments
samtools view -h pb-shuf.bam > pacbio_sup_sec.sam

minimap2 -ax sr --MD ../references/toy_reference.fa ../references/toy_reads.fa | samtools sort -o toy.bam -
samtools index toy.bam

minimap2 -ax map-hifi --MD --end-bonus 5  $AAV2 ../reads/aav2_subs.fa | samtools sort -o aav2_subs.bam - 
samtools index aav2_subs.bam

minimap2 -ax map-hifi --MD $AAV2 -B 1.5 --end-bonus 5 ../references/AAV2_AAV3.fa | samtools sort -o aav23.bam -
samtools index aav23.bam

minimap2 -ax map-hifi --MD $AAV2 -B 1.5 --end-bonus 5 ../references/wt2n496d389dna.fasta | samtools sort -o aav2389.bam -
samtools index aav2389.bam