# A quick-start guide to Nanopore assembly
This guide is intended as a brief overview guiding you through steps of generating an assembly from single molecule real-time (SMRT) sequencing using Oxford Nanopore, from initial sample preparation to a finished assembly. A more detailed summary and comparison of methods can be found at XXX.

## Pipeline overview:  
The general steps for an assembly using Oxford Nanopore Technology (ONT) is as follows:

Library prep --> Basecalling --> Assembly --> Assembly correction/polishing

There are example SLURM scripts provided below for each step of the process that can be run on the Cannon computing cluster at Harvard, or any equivalent cluster with a SLURM environment. All the software mentioned in this guide can be downloaded individually, ideally using `Conda` for ease of installation. Alternatively, a `Snakemake` pipeline that handles assembly and outputs a final, polished assembly is forthcoming :)

### Library Preparation  
Perhaps the most important step for generating high-quality, long read length sequences is obtaining high molecular weight, undamaged DNA from your samples. Using protocols (e.g. phenol-chloroform extraction) that avoid shearing/nicking your samples can result in read length N50 values > 20kb! Sample input requirements are higher compared to second generation technologies such as Illumina: about 1ug of DNA is recommended. Example protocols are available from the Bauer sequencing core.

Another important consideration is depth of coverage. To generate assemblies built from only ONT data, **~20-30X coverage** is sufficient as a rule of thumb for most assembly software. Hybrid assemblies built from both Illumina reads and ONT can use lower coverage: even **~5X coverage** ONT reads can greatly improve assembly contiguity, though a different assembler (e.g. MaSuRCA) will be needed; contact the Informatics department for details.

### Basecalling  
Output from ONT sequencers is in the form of FAST5 files, which is raw electrical signal from bases passing through pores on the membrane. There are numerous basecallers developed to translate raw signal to FASTQ format. For general purposes the Guppy basecaller, the "official" basecaller which is currently under active development by ONT, is the fastest, accurate, and is relatively easy to run.

For the latest ONT sequencers (e.g. the PromethION available at Harvard's sequencing core), the reads are basecalled by the sequencer itself, *so there is no need to run Guppy*. However, if you are using a different ONT sequencer (e.g. MinION) or if you wish to re-call the reads using an alternate model, you can download a binary executable for Guppy can from ONT's website. An account with the company is required for downloads, however. The GPU version is recommended, as it is much faster than the CPU version (e.g. *Phlox* reads basecalled in < 1 day, vs ~1 week!). If you do not have access to the download site, contact the Informatics department. An example script for basecalling using the GPU version of Guppy on the cluster is as follows:
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
Picking which assembler to use on your base-called data can be daunting due to the large number of options! In addition, you choice will depend on your organism and your needs for your assembly, as there is no "one size fits all" assembler: how big is your genome? Do you need as high accuracy as possible, or is the general scaffolding sufficient? In our experience, the `Flye` assembler generates contiguous assemblies with good base accuracy (after polishing steps, see below) and is far faster than any of the other assemblers we tested, making it a good place to start. Part of its speed is because `Flye` does not pre-correct the reads before assembly (compared to other assemblers e.g. `Canu`), which makes post-assembly polishing essential. Here is an example SLURM script for `Flye` using basic parameters:

<pre><code>#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o flye_dmel_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e flye_dmel_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python

#I had set up a Conda environment that I installed flye in. Change name as needed
source activate flye

#Replace with path to basecalled reads in FASTQ format
reads="/n/boslfs02/LABS/informatics/dkhost/diso1_a2_pass.fastq.gz"

#--nano-raw: indicates uncorrected Nanopore reads. Options for corrected reads.
#-g: Approximate genome size.
#-o: path where assembly will be output
#--threads: number of threads used for assembly
flye --nano-raw $reads -g 180m -o /n/boslfs02/LABS/informatics/dkhost/dmel_flye --threads 12

</code></pre>

A more detailed comparison of different assemblers, as well as a description of the their methods can be found at XXX. Scripts for running other assemblers are available at XXX.

### Polishing
Finally, due to the high error rate of ONT reads, assemblies need to be polished in order to achieve accuracies comparable to 2nd generation sequencing such as Illumina. This can be done in two ways: self-correction with ONT data (assuming high enough coverage) or correction with Illumina reads. Doing both polishing steps is ideal, though in our experience polishing with Illumina is far more impactful in improving accuracies. In addition, we have found that correction with Illumina *followed by* self-correction with ONT can lead to a **decrease** in accuracy, so if performing both rounds of polish it is important to self-correct *first*.

Self-correction is done by mapping the raw ONT reads to the assembly using programs such as `minimap2`, then calling a consensus. Several programs exist, however our current recommendation is to use `medaka` (which is quicker more accurate than other programs such as `nanopolish` and `racon`). It uses `minimap2` as a part of its run to map reads to the assembly, so there is no need to map them independently. Here is an example mapping to the *Drosophila melanogaster* assembly created above:

<pre><code>#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o medaka_dmel_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e medaka_dmel_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python

#I had set up a Conda environment that I installed medaka in. Change name as needed
source activate medaka

#Path to the raw ONT reads
reads="/n/boslfs02/LABS/informatics/dkhost/diso1_a2_pass.fastq"

#-i: raw read input
#-d: path to initial, unpolished assembly
#-o: directory to output to. Will be created if does not exist
#-t: number of threads
medaka_consensus -i $reads -d ./assembly.fasta -o dmel_flye.medaka -t 8
</code></pre>

Similar to self-correction, correction with Illumina involves first mapping the reads (ideally > ~40X coverage) to the assembly, which can be done with `minimap2` (or `bwa`/`bowtie2`), then calling a consensus. Consensus-calling can be done with several tools (e.g. `racon`, `pilon`), but our current recommendation is a recently developed tool known as `hypo` due to its increased accuracy and greatly increased speed relative to older tools. Here is an example using the *Drosophila melanogaster* assembly that we self-corrected using `medaka` above:

<pre><code>#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o hypo_dmel.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e hypo_dmel.err  # File to which STDERR will be written, %j inserts jobid

module load python samtools minimap2/2.9-fasrc01

#Path to the self-corrected assembly
ref="/n/boslfs02/LABS/informatics/dkhost/dmel_flye/dmel_flye.medaka/consensus.fasta"

#Path to the Illumina reads in FASTQ format
read1="/n/boslfs02/LABS/informatics/dkhost/SRR6702604_1.fastq.gz"
read2="/n/boslfs02/LABS/informatics/dkhost/SRR6702604_2.fastq.gz"

#Mapping the Illumina reads using Minimap2
minimap2 -ax sr -t 12 $ref $read1 $read2 > dmel_SRR6702604_vs_flye.medaka.sam

#Filtering, sorting, and indexing the mapping file using samtools
samtools view -h -Sb -F 4 -o dmel_SRR6702604_vs_flye.medaka.mapped.bam dmel_SRR6702604_vs_flye.medaka.sam
samtools sort dmel_SRR6702604_vs_flye.medaka.mapped.bam dmel_SRR6702604_vs_flye.medaka.mapped.sorted
samtools index dmel_SRR6702604_vs_flye.medaka.mapped.sorted.bam

#I had set up a Conda environment that I installed hypo in. Change name as needed
source activate hypo

#Hypo requires the Illumina read paths in a text file, one file per line
echo -e "$read1\n$read2" > il_names.txt

#-d: reference genome to polish
#-r: file with paths to Illumina reads
#-s: approximate genome size
#-c: approximate mean coverage of the Illumina reads over the assembly
#-b: sorted BAM file of Illumina reads against assembly
#-t: number of threads to use
#-o: output file to create
hypo -d $ref -r @il_names.txt -s 160m -c 30 -b dmel_SRR6702604_vs_flye.medaka.mapped.sorted.bam -t 12 -o dmel_flye.medaka.hypo.fasta
</code></pre>

It is possible to do multiple rounds of polishing; however, the first round is the most important, after which there are diminishing returns.

Using 30X ONT data coupled with ~30X Illumina data, assembly with `Flye` plus a single round of polishing with `hypo` generated an assembly with a BUSCO completeness score of **98.6%** (using the *Dipteran* gene data set) that was highly contiguous: almost entire chromosome arms were contained in single scaffolds. We also compared methods using a slightly larger, more repeat-dense insect genome, *Nasonia vitripennis*, where the improvements of `Flye` were even more noticeable: N50 of `Flye` assembly was an order of magnitude better than other assemblies!
