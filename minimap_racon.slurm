

module load python

source activate longread_toolkit

target="/scratch/dkhost/tetra_wtdbg2/tetra_wtdbg2_assem.ctg.fa"
reads="/scratch/lmcai/00_nanopore_raw_fastq/Tetrastigma_nanopore.fastq"

minimap2 -x map-ont -t 8 $target $reads | gzip -1 > tetra_wtdbg2_assem.ctg.minimap.paf.gz

racon -t 8 $reads tetra_wtdbg2_assem.ctg.minimap.paf.gz $target
