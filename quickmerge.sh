

prefix="tetra_canu_vs_miniasm_qm"

ref="/n/scratchlfs/informatics/dkhost/tetra.miniasm.pilon2_cleaned.fasta"
query="/n/scratchlfs/informatics/dkhost/tetrastigma.contigs_cleaned.fasta"

opt_l="495000" #set value for -l argument
#-l: controls length cutoff for anchor contigs. N50 is rule of thumb
#-ml: minimum alignment length for merging
options="-hco 5.0 -c 1.5 -l $opt_l -ml 20000"

~/quickmerge/MUMmer3.23/nucmer -l 100 -prefix $prefix $ref $query
echo "FINISHED NUCMER"
~/quickmerge/MUMmer3.23/delta-filter -i 90 -l 10000 -r -q $prefix.delta > $prefix.rq.delta
echo "FINISHED FILTER"
quickmerge -d tetra_canu_vs_miniasm_qm.rq.delta -q /n/scratchlfs/informatics/dkhost/tetrastigma.contigs_cleaned.fasta -r /n/scratchlfs/informatics/dkhost/tetra.miniasm.pilon2_cleaned.fasta -hco 5.0 -c 1.5 -l 495000 -ml 20000 -p $prefix
echo "FINISHED QUICKMERGE"

#ALTERNATIVELY:
#merge_wrapper.py -pre $prefix -hco 5.0 -c 1.5 -l $opt_l -ml 20000 $query $ref
