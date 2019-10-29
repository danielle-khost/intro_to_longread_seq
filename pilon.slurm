
module load python bowtie2 samtools jdk/1.8.0_172-fasrc01

#source activate assembly_eval

#readdir="/scratch/lmcai/00_tetrastigma_rafflesia_illumina_fastq/"

#r1="V5_S1_L001_R1_001"
#r2="V5_S1_L001_R2_001"
#r3="V5_S1_L002_R1_001"
#r4="V5_S1_L002_R2_001"

index="/scratch/dkhost/tetra_merge_v2/bowtie2_tetra_merge.racon/tetra_merge_v2.racon"
reads="/scratch/lmcai/20_tetrastigma_denovo_assembly/minimap_miniasm_pilon/tetra.illumina.fq.gz"
genome="/scratch/dkhost/tetra_merge_v2/merged_tetra_merge_v2.racon.fasta"

cd /scratch/dkhost/tetra_merge_v2/bowtie2_tetra_merge.racon

bowtie2-build $genome tetra_merge_v2.racon &> ind.out

cd /scratch/dkhost/tetra_merge_v2/pilon

#trim_galore --illumina --gzip --paired $readdir/${r1}.fastq.gz $readdir/${r2}.fastq.gz $readdir/${r3}.fastq.gz $readdir/${r4}.fastq.gz

bowtie2 -p 12 -x $index -U $reads | samtools view -o V5_vs_tetra_merge_v2.racon.mapped.bam -bS -h -F 4 -

samtools sort -o V5_vs_tetra_merge_v2.racon.mapped.sorted.bam V5_vs_tetra_merge_v2.racon.mapped.bam
samtools index V5_vs_tetra_merge_v2.racon.mapped.sorted.bam

pilon -Xmx1000G --genome $genome --unpaired V5_vs_tetra_merge_v2.racon.mapped.sorted.bam --output merged_tetra_merge_v2.racon.pilon --diploid --vcf --tracks --changes --threads 2
