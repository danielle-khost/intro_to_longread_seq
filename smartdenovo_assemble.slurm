#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-24:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o wtdbg2_ecoli1_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e wtdbg2_ecoli1_%j.err  # File to which STDERR will be written, %j inserts jobid

smartdenovo.pl -p ecoli_VREC0693.smartdenovo.assem.1 -c 1 ecoli_VREC0693.fastq.gz > ecoli_VREC0693.smartdenovo.assem.1.mak
make -f ecoli_VREC0693.smartdenovo.assem.1.mak
