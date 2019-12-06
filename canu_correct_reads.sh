#module load python

#source activate longread_toolkit

canu -correct -p dmel_corr_reads -d dmel_canu_corr genomeSize=140m useGrid=false -nanopore-raw /scratch/dkhost/diso1_a2_pass.fastq.gz
