# Long-read sequencing primer
## Introduction
Despite advances in sequencing technology, second-generation sequencing using short reads (such as Illumina or 454) still faces difficulties in certain areas of the genome such as regions with high repeat content, structural variants (e.g. CNVs), or sequences prone to PCR bias (e.g. high GC content). Recently developed so-called third-generation sequencing techniques attempt to address these shortcomings. The two dominant companies currently pioneering these techniques are Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT), both of which develop long-read, single molecule real-time sequencing. One limitation of these methods is the high error rate in the reads (~15% for Nanopore, ~13-15% for PacBio using latest chemistries), though this has been steadily improving and different methods have been developed to deal with error-prone reads (see below). While the two technologies share some broad stroke similarities, their methods and applications differ.

## How does it work?
### ONT MinION
A constant electrical current is run over proteinaceous membrane with tiny pores in it, creating an ion current. High molecular weight DNA is added to the membrane and the strands pass thru the pores, disrupting the current in a different manner depending on which bases are moving through the pore. Changes in the current are translated into nucleotide sequence in the base-calling step. Recently ONT released the MinION sequencer, which is small enough to be plugged into a standard USB port, and retails for ~$1000 USD.

#### Reads and library prep  
Double-stranded DNA is unzipped by a motor protein and a single strand is passed through the pore, producing a “1D” read with ~85% read accuracy in latest chemistries. Recently ONT has developed 1D2 reads: both the template strand and its reverse strand are sequenced one after the other without being physically ligated to each other, which they claim increases read accuracy to ~97%, though not all molecules can be sequenced in this manner.
•	MinKNOW software produces one file per DNA molecule in FAST5 format
•	Read length NG50 values range from 1-100kb. Ultra-long reads > 1Mb have also been reported

#### Base-calling  
Translates raw electrical signal of bases passing through a pore into nucleotide sequence. With current chemistries (R9.4), ~5 bases are passing through at a given time, leading to a large number of possible combinations of bases, coupled with inherently noisy signal.
-	Modern base-callers using trained neural networks
-	Choice of base-caller can have large effect on read accuracy (Wick 2019, Rang 2018)
-	Training base-caller on data similar to your organism of interest can improve accuracy (Wick 2019)
-	Current base-callers:
-	Developed by ONT:
-	Albacore: general purpose base-caller, available to ONT customers. Improved accuracy by adding transducer and reading raw electrical signal directly.
-	Guppy: faster version of Albacore, currently under active development. Has alternate models for base-calling, and can be trained with custom models.
-	Scrappie and Flappie: open-source base-callers. Flappie being actively developed.
-	3rd party:
-	Chiron

### PacBio
PacBio sequencing is slightly more developed than ONT, with earliest chemistries being introduced around 2009. Film with tiny cavities ten of nanometers in volume (“zero mode waveguides”, ZMGs) is affixed to glass plate. Each ZMW has a DNA polymerase anchored to the bottom of it, along with a single circular DNA molecule with insert, adapter and barcodes. A laser shone through the plate excites fluorescent bases incorporated by the polymerase as it synthesizes a new strand. PacBio’s most recent chemistry is called Sequel, with the previous generation called RS II.

-	Sequencing starts at the adapter and goes in a loop, potentially reading over the insert twice or more
-	Groups have reported that highly accurate reads (>99%) can be generated from taking consensus of circular subreads (Wenger et al 2019)
-	Risk of a ZMW containing no template (all background) or multiple (intercalated reads), making filtering necessary to remove low-quality reads
-	Read lengths vary from 1-15kb, with reads > 90kb possible
-	Read accuracies around ~86%


## What can I do with it?
The most common applications of long read sequencing is to be able to directly study regions of the genome that are still problematic for 2nd generation technologies, i.e. highly repetitive regions such as the pericentromeric heterochromatin, multi-copy gene families such as MHC genes, regions of structural variation, or sequences known to be affected by PCR bias. Previously the high cost of 3rd generation long read sequencing limited its application to organisms with small genomes. However, improvements in data output and error rate has allowed accessible sequencing of more complex genomes. Recent examples the kinds of questions you can answer using long reads:

-	Characterize structure of tandem repeat arrays (e.g. satellite DNA)
-	Ultra-long Nanopore reads used to describe properties of satellite arrays in grass pea plant species with large genome enriched for satellites (6 Gb) (Vondrak 2019)
-	PacBio used to characterize structure of complex satDNA loci in D. melanogaster (Khost 2017)
-	Detecting structural variation and copy number variation
-	Use circular consensus sequencing protocol (multiple sequencing passes over single template) for PacBio to generate extremely accurate long reads to identify SNPs, indels, and structural variants in human population data (Wenger et al 2019)
-	Designed algorithms to account for high error rate and align long reads to identify structural variants, which is able to identify SVs difficult to identify with short reads such as nested SVs, inverted tandem duplications, or inversions flanked by indels (Sedlazeck 2018)
-	Do comparative genomics in difficult to assemble regions
-	Sequenced several great ape genomes using high-coverage long reads plus Illumina sequencing and identified structural variants unique to each lineage, as well as compared the differences in repetitive element landscapes (Kronenberg 2018)

## What do I need?
Despite improvements, error rates of ONT and PacBio sequencing remain too high to assemble the reads straight away, necessitating an error-correction of the reads and/or assembly. One important consideration when designing a long-read sequencing experiment is whether you will be using only long reads for your de novo assembly (self-correction) or whether supplementary Illumina reads will be used as well (hybrid assembly).
-	If using self-correction (long reads only):
-	PacBio: to obtain accuracies close to Illumina, need very high coverage (>90X)
-	Nanopore: ~30X coverage is sufficient for self-correction. However, still has accuracy issues and biases in homopolymer regions (see below)
-	If doing hybrid assembly:
-	Even low coverage long reads (5-10X) can greatly improve the contiguity of short read assemblies
-	For correcting assemblies using Pilon

## What assemblers can I use?
- *Minimap/miniasm*: Read mapper and assembler, respectively. Specialized for fast assembly of noisy ONT reads, especially in larger genomes. Note: does not attempt to assemble repetitive regions! Also does not perform an initial error correction of reads like other long read assemblers, resulting in lower per-base accuracy.
- *SMARTdenovo*: Another rapid assembler usable for ONT and PacBio. Like miniasm, it does not perform error-correction, instead doing an all-by-all alignment of raw reads. Consensus polishing tools (e.g. Quiver for PacBio or Nanopolish for ONT) are required for accurate assembly. Can also perform error correction using other assemblers (e.g. Canu), then pass corrected sequences to SMARTdenovo.
- *Canu*: Assembler specialized for noisy long-read sequences, useable for PacBio or ONT. First finds overlaps between reads using a hashing algorithm (MHAP) for faster alignment, then generates a corrected consensus based on overlaps and performs an assembly using the corrected sequences. Due to the correction step, it is more computationally intensive that miniasm but is still fairly fast.
- *Falcon*: PacBio only. Diploid-aware, allowing phasing of genotypes with sufficient coverage. Specialized for larger genomes. Similar to Canu, performs an error correction step where smaller reads are aligned against largest subset of reads to generate corrected consensus reads that are used for assembly. In addition to primary contigs, also produces “haplotigs” in regions of divergent haplotypes.
- *Shasta*: ONT only, recently developed in response to heavy computational requirements for human ONT assembly. Developed to be fast (human assembly < 1 day) and more resilient to errors in homopolymer regions that ONT has trouble with. Loads every single base into memory all at once, which for large genomes can require a HUGE amount of processing power, e.g. 128 CPUs and ~1950 GB of RAM for a human genome at 60X coverage!!
- *MaSuRCa*: Useable for PacBio, ONT and Illumina data. Allows hybrid assembly of long reads: combines short Illumina reads together to form super-reads, which are then approximately aligned using k-mers (similar to MHAP algorithm) with long reads to produce long, accurate “mega-reads”. These reads are then used for assembly using a modified Celera assembler (CABOG). MaSuRCa has been used on difficult, repeat-dense genomes such as *A. tauschii*, a plant related to wheat.


## How do I correct the reads?
To compensate for the high error rates of the reads, need to either use self-correction or using short reads (e.g. from Illumina). Correction can occur before assembly on the reads themselves, on the finished assembly, or both.
- *Quiver/Arrow*: PacBio only, generates a consensus using self-correction of PacBio reads: the longest subset of reads is corrected by aligning shorter reads. Used after assembly is constructed as a final polishing step, can be done in conjunction with other polishing software (e.g. Pilon).
- *Pilon*: Correction using short reads, performed on finished assembly. Input a BAM file of short reads mapped against genome, which allows Pilon to correct SNPs and indels, as well as fill gaps and identify local mis-assemblies.
- *Racon*: PacBio and ONT. Intended for use with assemblers that do not perform an initial error correction step (e.g. miniasm, SMARTdenovo). Can be used for correction of the reads or as final assembly polishing.
- *Nanopolish*: self-correction specialized for ONT data, run on finished assembly. Unlike other polishing programs, requires the raw signal-level data from ONT in addition to the output from the base-calling software. Base-called reads are mapped against the draft assembly, which is input to Nanopolish.

## Should I use PacBio or Nanopore?
Both PacBio and ONT are very active being developed, and as such the pros and cons of each technology can change rapidly. Currently, there is no clear-cut answer as to which technology is “better,” and which one you choose will depend on your needs and resources. However, there are several things to keep in mind:
-	Cost: PacBio sequencers are prohibitively expensive for most labs, making it necessary to rely on sequencing cores. ONT’s MinION sequencer, however, can be purchased for $1000, and produces a similar amount of data per flow cell.
-	Genome size: while long read sequencing is improving data output and creating more efficient assembly algorithms, and there have been several assemblies of complex genomes published, many protocols are still geared towards smaller genomes.
  -	For smaller genomes (e.g. bacteria, some invertebrates, etc.), either technology is fairly accessible and can be assembled using standard benchtop computers. ONT’s MinION is probably the best route due to the high output from a single flow cell and low cost.
  -	For larger genomes, likely need access to computing cluster to assemble in practical amount of time
-	ONT offers higher data output, but groups have reported error rate is still a problem and PacBio produces better results
-	Accuracy:
  -	One recent study found that even with custom trained base-calling model and polishing with Nanopolish, the max accuracy was 99.94% for ONT
  -	Unlike PacBio, ONT errors are not randomly distributed: it has trouble with long homopolymers due to uncertainty in base-calling. While recent base callers have improved this bias, some groups still report it as an issue.
  -	Most groups, even when using high coverage sequencing, still correct their assemblies with Illumina
-	Material: PacBio require much more input for sequencing
-	High error rate of MinION is fine for SV detection, but if looking at SNPs or indels or haplotype variation, accuracy is still problem
-	**Take-away: ONT produces longer reads and more data for less money, at the cost of an increased error rate; PacBio is more expensive and produces shorter reads, but has a better accuracy**

## Example pipelines
Which assembly process you use will depend on your organism and a variety of other factors, and each method has a wide variety of advanced usage. However, here are some basic protocols for assembly ONT reads using example data from *E. coli* to use as a foundation:  
- *canu_ecoli_assemble.slurm*: performs full Canu pipeline (correction, trimming, assembly) using default parameters. Disables grid usage (useGrid=false), restricting resources to single local machine (e.g. a personal computer or lab server)...this should be fine for small genomes (like bacteria), but larger genomes should use the grid.
- *canu_ecoli_assemble_grid.slurm*: as above, but will run on the grid (i.e. Odyssey). By default, Canu will auto-detect available resources and will configure itself; however, it does not request explicit time limits or partitions. You can specify these using the gridOptions=" " argument, and can also control resources used by individual steps (e.g. correction, overlapping, etc). **WIP**
- *miniasm_ecoli_assemble.slurm*: performs overlapping and assembly using minimap2 and miniasm, plus self-correction of the assembly using Racon. This should run quickly, but is error-prone: further correction with Illumina reads, or pre-correction of ONT reads is useful. **WIP**
