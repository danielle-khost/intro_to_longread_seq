#!/bin/bash

#Explanation of arguments:
#-p: prefix that will be appended to output files
#-d: creates directory that will store intermediate and output files
#genomeSize: rough estimate of total genome size of the organism
#-nanopore-raw: specifies that reads are uncorrected nanopore reads. Can also give corrected reads
#ecoli_VREC0693.fastq.gz: FASTQ file from basecalling step containing raw reads. Can also be in FASTA format

module load canu

canu -p ecoli -d ecoli-ONT_grid2 genomeSize=4.8m useGrid=true gridOptions="--time=24:00:00 --partition=general" -nanopore-raw ecoli_VREC0693.fastq.gz
