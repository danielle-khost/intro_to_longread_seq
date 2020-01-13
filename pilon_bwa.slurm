
module load python samtools jdk/1.8.0_172-fasrc01 bwa

#source activate assembly_eval

readsF="/n/holylfs/LABS/srivastava_lab/genome_raw_reads_FINAL/Reddien_lanes_L3_L4_L5_combined_F_FINAL.fastq"
readsR="/n/holylfs/LABS/srivastava_lab/genome_raw_reads_FINAL/Reddien_lanes_L3_L4_L5_combined_R_FINAL.fastq"
genome="/scratch/dkhost/hofsenia_shasta/ShastaRun/hofsenia_shasta.fasta"

cd /scratch/dkhost/hofsenia_shasta/ShastaRun/
echo "START"
date
bwa index -p hofsenia_shasta $genome
echo "INDEXED GENOME"
date

cd /scratch/dkhost/hofsenia_shasta/pilon/

echo "START MAPPING"
date
bwa mem -t 18 $genome $readsF $readsR > reddien_vs_hofsenia_shasta.sam
echo "MAPPED READS"
date

samtools view -o reddien_vs_hofsenia_shasta.mapped.bam -bS -h -F 4 reddien_vs_hofsenia_shasta.sam
samtools sort -o reddien_vs_hofsenia_shasta.mapped.sorted.bam reddien_vs_hofsenia_shasta.mapped.bam
samtools index reddien_vs_hofsenia_shasta.mapped.sorted.bam
echo "SORTED AND INDEXED"

echo "START PILON"
date
pilon -Xmx900G --genome $genome --frags reddien_vs_hofsenia_shasta.mapped.sorted.bam --output hofsenia_shasta.pilon --diploid --vcf --tracks --changes --threads 2
echo "END PILON"
date
