# Background on long-read sequencing
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
-	Guppy: faster version of Albacore, currently under active development. Has alternate models for base-calling, and can be trained with custom models.
-	Albacore: general purpose base-caller, available to ONT customers. Improved accuracy by adding transducer and reading raw electrical signal directly.
-	Scrappie and Flappie: open-source base-callers. Flappie being actively developed.

There are also 3rd party basecallers (e.g. Chiron) being developed, though they are less commonly used.

### PacBio
PacBio sequencing is slightly more developed than ONT, with earliest chemistries being introduced around 2009. Film with tiny cavities ten of nanometers in volume (“zero mode waveguides”, ZMGs) is affixed to glass plate. Each ZMW has a DNA polymerase anchored to the bottom of it, along with a single circular DNA molecule with insert, adapter and barcodes. A laser shone through the plate excites fluorescent bases incorporated by the polymerase as it synthesizes a new strand. PacBio’s most recent chemistry is called Sequel, with the previous generation called RS II.

-	Sequencing starts at the adapter and goes in a loop, potentially reading over the insert twice or more
-	Groups have reported that highly accurate reads (>99%) can be generated from taking consensus of circular subreads (Wenger et al 2019)
-	Risk of a ZMW containing no template (all background) or multiple (intercalated reads), making filtering necessary to remove low-quality reads
-	Read lengths vary from 1-15kb, with reads > 90kb possible
-	Read accuracies around ~86%
