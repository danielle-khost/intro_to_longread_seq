#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-12:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p general   # Partition to submit to
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o miniasm_ecoli_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e miniasm_ecoli_%j.err  # File to which STDERR will be written, %j inserts jobid

module load minimap2 miniasm racon

#Generate all-by-all overlap of reads
minimap2 -x ava-ont -t 8 ecoli_VREC0693.fastq.gz ecoli_VREC0693.fastq.gz | gzip -1 > ecoli_VREC0693.paf.gz

#Generate assembly
miniasm -f ecoli_VREC0693.fastq.gz ecoli_VREC0693.paf.gz > ecoli_VREC0693.assem.gfa

#Correction with racon
##Convert from GFA to FASTA
awk '$1 ~/S/ {print ">"$2"\n"$3}' ecoli_VREC0693.assem.gfa > ecoli_VREC0693.assem.fasta

##Map reads against alignment
minimap2 -t 8 ecoli_VREC0693.assem.fasta ecoli_VREC0693.fastq.gz > ecoli_VREC0693.assem.gfa.paf

##Take consensus using racon
racon -t 8 ecoli_VREC0693.fastq.gz ecoli_VREC0693.assem.gfa.paf ecoli_VREC0693.assem.fasta ecoli_VREC0693.assem.racon.fasta
