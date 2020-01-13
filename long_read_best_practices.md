# Best practices guide to Nanopore assembly
This guide is intended as a brief overview guiding you through steps of generating an assembly from single molecule real-time (SMRT) sequencing using Oxford Nanopore, from initial sample preparation to a finished assembly. A more detailed summary and comparison of methods can be found at XXX.

All the software mentioned in this guide can be downloaded individually. Alternatively, the latest versions are packaged in a Conda environment here (link_to_env) for easy installation and running on the Cannon computing cluster.

## Pipeline overview:  
The general steps for an assembly using Oxford Nanopore Technology (ONT) is as follows:

Library prep --> Basecalling --> Assembly --> Assembly correction/polishing

### Library Preparation  
Perhaps the most important step for generating high-quality, long read length sequences is obtaining high molecular weight, undamaged DNA from your samples. Using protocols (e.g. phenol-chloroform extraction) that avoid shearing/nicking your samples can result in read length N50 values > 20kb! Sample input requirements are higher compared to second generation technologies such as Illumina: about 1ug of DNA is recommended. Example protocols are available from the Bauer sequencing core.

Another important consideration is depth of coverage. To generate assemblies built from only ONT data, **~20-30X coverage** is sufficient as a rule of thumb for most assembly software. Hybrid assemblies built from both Illumina reads and ONT can use lower coverage: even **~5X coverage** ONT reads can greatly improve assembly contiguity.

### Basecalling  
Output from most ONT sequencers (e.g. MinION) is in the form of FAST5 files, which is raw electrical signal from bases passing through pores on the membrane. There are numerous basecallers developed to translate raw signal to FASTQ format. For general purposes the Guppy basecaller, the "official" basecaller which is currently under active development by ONT, is the fastest, accurate, and is relatively easy to run.

A binary executable for Guppy can be downloaded from ONT's website, though an account with the company is required. The GPU version is recommended, as it is much faster than the CPU version (e.g. *Flox* reads basecalled in < 1 day, vs ~1 week!). If you do not have access to the download site, contact the Informatics department. An example script for basecalling using the GPU version of Guppy on the cluster is as follows:
<pre><code>#!/bin/bash

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

fast5="/path/to/fast5/files"
outdir="/path/to/outdir"

~/ont-guppy/bin/guppy_basecaller -i $fast5 \ #Change the path to where you installed the executable
 -s $outdir \
--flowcell FLO-MIN106 \ #This value depends on which flow cell was used to generate the data, in this case a MinION r9.4 flow cell
--kit SQK-LSK109 \ #This value depends on the kit used
--compress_fastq \
--device cuda:0

</code></pre>

### Assembly
Picking which assembler to use on your base-called data can be daunting due to the large number of options! In addition, you choice will depend on your organism and your needs for your assembly, as there is no "one size fits all" assembler: how big is your genome? Do you need as high accuracy as possible, or is the general scaffolding sufficient? In our experience, the `wtdbg2` assembler generates contiguous assemblies with good base accuracy (after polishing steps, see below) while still be computationally efficient enough to assemble large (>1gb) genomes in under a few days. Part of its speed is because `wtdbg2` does not pre-correct the reads before assembly (compared to other assemblers e.g. `Canu`), which makes post-assembly polishing essential. Here is an example SLURM script:

<pre><code>#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-48:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o wtdbg2_dmel1_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e wtdbg2_dmel1_%j.err  # File to which STDERR will be written, %j inserts jobid

#Assembler step
#-x: uses nanopore presets
#-g approx. genome size
#-i: input raw basecalled sequence
#-t: number of threads
#-fo: force output to files with following prefix
wtdbg2 -x nanopore -g 140m -i ../diso1_a2_pass.fastq.gz -t 8 -fo dmel_wtdbg2.assem.1

#Consensus calling step
wtpoa-cns -t 8 -i dmel_wtdbg2.assem.1.ctg.lay.gz -fo dmel_wtdbg2.assem.1.ctg.fa

</code></pre>

A more detailed comparison of different assemblers, as well as a description of the their methods can be found at XXX. Scripts for running other assemblers are available at XXX.

### Polishing
Finally, due to the high error rate of ONT reads, assemblies need to be polished in order to achieve accuracies comparable to 2nd generation sequencing such as Illumina. This can be done in two ways: self-correction with ONT data (assuming high enough coverage) or correction with Illumina reads. Doing both polishing steps is ideal, though polishing with Illumina is far more impactful in improving accuracies. In addition, we have found that correction with Illumina *followed by* self-correction with ONT can lead to a **decrease** in accuracy, so if performing both rounds of polish it is important to self-correct *first*. Self-correction is done by first mapping the ONT reads to the assembly using `minimap2`, then polishing with programs such as `nanopolish` or `racon`. An example script using `minimap2` and `racon` (which is quicker than `nanopolish`) is as follows:

<pre><code>target="/scratch/dkhost/dmel_wtdbg2/dmel_wtdbg2.pilon.fasta"
reads="/scratch/dkhost/diso1_a2_pass.fastq.gz"

minimap2 -x map-ont -t 8 $target $reads | gzip -1 > dmel_wtdbg2.pilon.racon.paf.gz

~/racon/build/bin/racon -t 8 $reads dmel_wtdbg2.pilon.racon.paf.gz $target
</code></pre>

Similar to self-correction, correction with Illumina involves first mapping the reads (ideally > ~40X coverage) to the assembly, which can be done with `minimap2` or also `bwa`/`samtools`, then consensus-calling using `pilon`:

<pre><code>readsF="/scratch/dkhost/SRR6702604_1.fastq.gz"
readsR="/scratch/dkhost/SRR6702604_2.fastq.gz"
genome="/scratch/dkhost/dmel_wtdbg2/dmel_wtdbg2.assem.1.ctg.fa"

cd /scratch/dkhost/dmel_wtdbg2/
echo "START"
date
bwa index $genome
echo "INDEXED GENOME"
date

echo "START MAPPING"
date
prefix="SRR6702604_vs_wtdbg2"
bwa mem -t 12 $genome $readsF $readsR > ${prefix}.sam
echo "MAPPED READS"
date

samtools view -o ${prefix}.mapped.bam -bS -h -F 4 ${prefix}.sam
samtools sort -o ${prefix}.mapped.sorted.bam ${prefix}.mapped.bam
samtools index ${prefix}.mapped.sorted.bam
echo "SORTED AND INDEXED"

echo "START PILON"
date
pilon -Xmx200G --genome $genome --frags ${prefix}.mapped.sorted.bam --output dmel_wtdbg2.pilon --diploid --vcf --threads 2
echo "END PILON"
date
</code></pre>

It is possible to do multiple rounds of polishing with `pilon`; however, the first round is the most important, after which there are diminishing returns.

Using 30X ONT data coupled with ~60X Illumina data, assembly with `wtdbg2` plus a single round of polishing with `pilon` generated an assembly with a BUSCO completeness score of **98.1%** (using the *Dipteran* gene data set) that was highly contiguous: almost entire chromosome arms were contained in single scaffolds.
