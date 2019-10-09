#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 6-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o qm_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e qm_%j.err  # File to which STDERR will be written, %j inserts jobid

prefix="tetra_canu_vs_miniasm_qm"

ref="/n/scratchlfs/informatics/dkhost/tetra.miniasm.pilon2_cleaned.fasta"
query="/n/scratchlfs/informatics/dkhost/tetrastigma.contigs_cleaned.fasta"

opt_l="495000" #set value for -l argument
#-l: controls length cutoff for anchor contigs. N50 is rule of thumb
#-ml: minimum alignment length for merging
options="-hco 5.0 -c 1.5 -l $opt_l -ml 20000"

~/quickmerge/MUMmer3.23/nucmer -l 100 -prefix $prefix $ref $query
~/quickmerge/MUMmer3.23/delta-filter -i 90 -l 10000 -r -q $prefix.delta > $prefix.rq.delta
~/quickmerge/quickmerge -d $prefix.rq.delta -q $query -r $ref $options

#ALTERNATIVELY:
#merge_wrapper.py -pre $prefix -hco 5.0 -c 1.5 -l $opt_l -ml 20000 $query $ref
