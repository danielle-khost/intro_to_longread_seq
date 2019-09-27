#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-24:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p general   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o canu_ecoli_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e canu_ecoli_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python

source activate longread_toolkit

./assemble.sh
