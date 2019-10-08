#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 6-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p general   # Partition to submit to
#SBATCH --mem=48G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o shasta_dmel_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e shasta_dmel_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python

source activate assembly_eval

run_BUSCO.py -f -n 12 -m geno -i /n/scratchlfs/informatics/dkhost/tetra_shasta/Assembly.fasta -o busco_tetra_shasta -l /n/scratchlfs/informatics/dkhost/embryophyta_odb9/
