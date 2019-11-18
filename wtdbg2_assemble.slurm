#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-48:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o wtdbg2_dmel1_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e wtdbg2_dmel1_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python
source activate longread_toolkit

#Assembler step
#-x: uses nanopore presets
#-g approx. genome size
#-i: input raw basecalled sequence
#-t: number of threads
#-fo: force output to files with following prefix
wtdbg2 -x nanopore -g 140m -i ../diso1_a2_pass.fastq.gz -t 8 -fo dmel_wtdbg2.assem.1

#Consensus calling step
wtpoa-cns -t 8 -i dmel_wtdbg2.assem.1.ctg.lay.gz -fo dmel_wtdbg2.assem.1.ctg.fa
