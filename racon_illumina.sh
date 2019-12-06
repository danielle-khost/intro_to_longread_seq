module load python samtools

#source activate longread_toolkit

target="/scratch/dkhost/tetra_merge_v2/merged_tetra_merge_v2.racon.fasta"
reads="/scratch/lmcai/20_tetrastigma_denovo_assembly/minimap_miniasm_pilon/tetra.illumina.fq.gz"

#minimap2 -x map-ont -t 8 $target $reads | gzip -1 > tetra_wtdbg2_assem.ctg.minimap.paf.gz

samtools view -h -o /scratch/dkhost/tetra_merge_v2/pilon/V5_vs_tetra_merge_v2.racon.mapped.sorted.sam /scratch/dkhost/tetra_merge_v2/pilon/V5_vs_tetra_merge_v2.racon.mapped.sorted.bam

~/racon/build/bin/racon -t 24 $reads /scratch/dkhost/tetra_merge_v2/pilon/V5_vs_tetra_merge_v2.racon.mapped.sorted.sam $target
