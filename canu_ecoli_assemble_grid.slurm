#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-12:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p general   # Partition to submit to
#SBATCH --mem=8G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o canu_ecoli_grid_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e canu_ecoli_grid_%j.err  # File to which STDERR will be written, %j inserts jobid

module load canu

#Run this script after running basecalling algorithm on raw signal output from sequencing (using Guppy, Albacore, etc.)...should produce a FASTQ file

#This is the most basic way to run canu, it uses all the default parameters and runs the whole pipeline: read correction, trimming, and assembly

#Explanation of arguments:
#-p: prefix that will be appended to output files
#-d: creates directory that will store intermediate and output files
#genomeSize: rough estimate of total genome size of the organism
#-nanopore-raw: specifies that reads are uncorrected nanopore reads. Can also give corrected reads
#ecoli_VREC0693.fastq.gz: FASTQ file from basecalling step containing raw reads. Can also be in FASTA format

canu -p ecoli -d ecoli-ONT_grid2 genomeSize=4.8m useGrid=true gridOptions="--time=24:00:00" -nanopore-raw ecoli_VREC0693.fastq.gz
