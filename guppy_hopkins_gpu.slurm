#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH -t 4-00:00                # Runtime in D-HH:MM
#SBATCH --partition=gpu
#SBATCH --mem=40GB
#SBATCH --gres=gpu:1
#SBATCH -J guppy
#SBATCH -o outfile.%A.out # File to which STDOUT will be written
#SBATCH -e outfile.%A.err # File to which STDERR will be written
#SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=danielle_khost@fas.harvard.edu # Email to which notifications will be sent

fast5="/n/boslfs/INSTRUMENTS/illumina/Nanopore/Hopkins/20191025_Hopkins_MP_409_R_SRE/Hopkis_MP_409_R_SRE/20191025_1851_MN27108_FAL55029_daa1a5c7/fast5/"
outdir="/n/holylfs/LABS/informatics/dkhost/guppy/20191025_Hopkins_MP_409_R_SRE"

~/ont-guppy/bin/guppy_basecaller -i $fast5 \
 -s $outdir \
--flowcell FLO-MIN106 \
--kit SQK-LSK109 \
--compress_fastq \
--device cuda:0
