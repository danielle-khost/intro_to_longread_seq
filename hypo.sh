#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o hypo_dmel_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e hypo_dmel_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python samtools
source activate longread_toolkit

ref="/n/boslfs02/LABS/informatics/dkhost/dmel_flye/dmel_flye.medaka/consensus.fasta"
read1="/n/boslfs02/LABS/informatics/dkhost/SRR6702604_1.fastq.gz"
read2="/n/boslfs02/LABS/informatics/dkhost/SRR6702604_2.fastq.gz"

minimap2 -ax sr -t 12 $ref $read1 $read2 > dmel_SRR6702604_vs_flye.medaka.sam

samtools view -h -b -F 4 -o dmel_SRR6702604_vs_flye.medaka.mapped.bam dmel_SRR6702604_vs_flye.medaka.sam
samtools sort -o dmel_SRR6702604_vs_flye.medaka.mapped.sorted.bam dmel_SRR6702604_vs_flye.medaka.mapped.bam
samtools index dmel_SRR6702604_vs_flye.medaka.mapped.sorted.bam

source activate hypo

echo -e "$read1\n$read2" > il_names.txt

hypo -d $ref -r @il_names.txt -s 160m -c 30 -b dmel_SRR6702604_vs_flye.medaka.mapped.sorted.bam -t 12 -o dmel_flye.medaka.hypo.fasta
