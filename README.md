# De novo assembly of Illumina reads using ABySS and alignment using BWA

This workshop is based on the original design by [Shaun Jackman](http://sjackman.ca).

# Purpose

We will use ABySS to assemble a 200 kbp region using paired-end reads from the Illumina platform. BWA-MEM is used to align the assembled contigs to the human reference genome. IGV is used to visualize these alignments and variants.

After this lab, you will have learned how to use ABySS to assemble a small genome, use BWA-MEM to align reads and contigs to a reference genome, and use IGV to visualize these alignments.

# Contents

* [Getting Started](#getting-started)
* [Exercise 0: Index the reference using BWA](#exercise-0-index-the-reference-using-bwa)
* [Exercise 1: Align the reads to the reference using BWA (optional)](#exercise-1-align-the-reads-to-the-reference-using-bwa-optional)
* [Exercise 2: Inspect the reads](#exercise-2-inspect-the-reads)
* [Exercise 3: Assemble the reads into contigs using ABySS](#exercise-3-assemble-the-reads-into-contigs-using-abyss)
* [Exercise 4: Align the contigs to the reference using web BLAT](#exercise-4-align-the-contigs-to-the-reference-using-web-blat)
* [Exercise 5: Align the contigs to the reference using BWA-MEM](#exercise-5-align-the-contigs-to-the-reference-using-bwa-mem)
* [Exercise 6: Browse the contig to reference alignments using samtools tview](#exercise-6-browse-the-contig-to-reference-alignments-using-samtools-tview)
* [Exercise 7: Browse the contig to reference alignments using IGV](#exercise-7-browse-the-contig-to-reference-alignments-using-igv)

# Getting Started

## System requirements

+ ABySS: genome sequence assembler for short reads
  http://www.bcgsc.ca/platform/bioinfo/software/abyss
+ BWA: align short reads and long contigs
  http://bio-bwa.sourceforge.net
+ IGV: visualize a genome
  http://www.broadinstitute.org/igv
+ samtools: manipulate SAM/BAM alignment files
  http://www.htslib.org

## Install the software using conda

Install ABySS, BWA, IGV, and samtools using conda.

```sh
conda install abyss bwa igv samtools
```

## Download the data

Create your working directory.

```sh
mkdir ~/wgs_tutorial
cd ~/wgs_tutorial
```

Download the FASTQ files.

```sh
curl -LO http://bioinfo.cnio.es/Cursos/Master2018/wgs/bacsim.read1.fastq.gz
curl -LO http://bioinfo.cnio.es/Cursos/Master2018/wgs/bacsim.read2.fastq.gz
```

Download the reference file of human chromosome 3.

```sh
curl -LO http://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr3.fa.gz
gunzip chr3.fa.gz
```

# Exercise 0: Index the reference using BWA

Index the reference file.

```sh
bwa index chr3.fa
```

~5 min

# Exercise 1: Align the reads to the reference using BWA (optional)

Run BWA.

```sh
bwa mem -t2 chr3.fa bacsim.read1.fastq.gz bacsim.read2.fastq.gz > bacsim.sam
```

~40 min

While this job is running in the background, open a new terminal and continue with the next section.

# Exercise 2: Inspect the reads

There are two FASTQ files for this set of paired-end Illumina data: one file for the forward reads and one file for the reverse reads. Look at the first few reads.

```sh
gunzip -c bacsim.read1.fastq.gz | head
```

How long are the reads? (hint: use `awk '{print length}'`)

> ```sh
> gunzip -c bacsim.read1.fastq.gz | awk '{print length}' | sed -n '2p'
> ```
> 50 bp

How many lines are there in total (both files)? (hint: use `wc -l`)

> ```sh
> gunzip -c bacsim.read1.fastq.gz bacsim.read2.fastq.gz | wc -l
> ```
> 40,019,248 lines

How many lines per read?

> 4 lines per read

How many reads are there in total (both files)?

> 40019248 / 4 = 10,004,812 reads

How many bases are sequenced?

> 10004812 * 50 = 500,240,600 bp

Assuming the BAC is 200 kbp, what is the depth of coverage?

> 500240600 / 200000 = ~2500 fold coverage

# Exercise 3: Assemble the reads into contigs using ABySS

Run the assembly (you may want to take a short break while it runs).

```sh
mkdir k48
ln -s ../bacsim.read1.fastq.gz ../bacsim.read2.fastq.gz k48/
abyss-pe -C k48 name=HS0674 k=48 s=200 v=-v in="bacsim.read1.fastq.gz bacsim.read2.fastq.gz" contigs 2>&1 | tee abyss.log
```

~15 min

Look at the option `-n,--dry-run` of abyss-pe. Its output is the
commands that ABySS will run for the assembly.

```sh
abyss-pe name=HS0674 k=48 s=200 in="bacsim.read1.fastq.gz bacsim.read2.fastq.gz" contigs -n
```

The assembly runs in three stages: assemble contigs without paired-end information, align the paired-end reads to the initial assembly, and merge contigs joined by paired-end information. You can instruct ABySS to stop after any of these stages. Use the `-n` option to see the commands for each stage.

```sh
abyss-pe name=HS0674 k=48 in="bacsim.read1.fastq.gz bacsim.read2.fastq.gz" unitigs -n
abyss-pe name=HS0674 k=48 in="bacsim.read1.fastq.gz bacsim.read2.fastq.gz" pe-sam -n
abyss-pe name=HS0674 k=48 in="bacsim.read1.fastq.gz bacsim.read2.fastq.gz" contigs -n
```

Once the assembly has completed, look at the assembled contigs using `less`. The option `-S` disables line wrapping.

```sh
less -S k48/HS0674-contigs.fa
```

How many contigs are longer than 100 bp?

> 5

What is the length of the longest contig (hint: use `awk '{print length}'`)?

> ```sh
> awk '{print length}' k48/HS0674-contigs.fa | sort -n | tail -n1
> ```
> 113,438 bp

What is the N50 of the assembly?

```sh
abyss-fac k48/HS0674-contigs.fa
```

> ```
> n   n:500 L50  min   N75    N50     N25     E-size  max     sum     name
> 15  4     1    5258  70187  113438  113438  89993   113438  199226  k48/HS0674-contigs.fa
> ```
> 113,438 bp

View the assembly log.

```sh
less -S abyss.log
```

What portion of the reads align to the assembly? (hint: search for "Mapped")

> Mapped 5025757 of 10004812 reads (50.2%)

# Exercise 4: Align the contigs to the reference using web BLAT

Open BLAT in a web browser: <http://genome.ucsc.edu>

View the assembled contigs in a text editor. `gview` is one text editor. You can use any text editor that you prefer.

```sh
gview k48/HS0674-contigs.fa
```

Disable line wrap to make it easier to select the full sequence. For Emacs, select the option "Options -> Line Wrapping in this Buffer -> Truncate Long Lines." For Vim, select the option "Edit -> File Settings -> Toggle Line Wrap", or type `:set nowrap`

Select the contig whose length is approximately 70 kbp and copy-and-paste its sequence into BLAT.

What is the exact length of this contig?

> 70,187 bp

Click "browser" for the best alignment and then zoom out 10x. To which chromosome and band does this contig align?

chr3:187,585,392-187,625,896

> chr3q27.3

What are the nearest two genes?

> SST and RTP2

Set the "Common SNPs" track in Variation to "pack".

A SNV between the assembled sequence and reference sequence is displayed with a red line. Zoom in on a SNV. Is it in dbSNP?

Set the "RepeatMasker" track in Repeats to "full". Look for "simple repeats" that overlap the contig. Did any of them cause problems for the assembly?

> Simple_repeat (TTCC)n [chr3:187597669-187597735]

Zoom in to see the sequence of the feature.

# Exercise 5: Align the contigs to the reference using BWA-MEM

Warning: do not start BWA-MEM until BWA has completed (if you decided to run it) unless your machine has at least 2 GB of RAM.

Run BWA-MEM.

```sh
bwa mem -t2 chr3.fa k48/HS0674-contigs.fa >k48/HS0674-contigs.sam
samtools sort -o k48/HS0674-contigs.bam k48/HS0674-contigs.sam
samtools index k48/HS0674-contigs.bam
```

~1 min

# Exercise 6: Browse the contig to reference alignments using samtools tview

Run samtools tview.

```sh
samtools tview k48/HS0674-contigs.bam chr3.fa
```

Go to the region `chr3:187,586,051`

```
g chr3:187,586,051
```

What variants do you see?

> a G/A SNV

Go to the region `chr3:187,597,641`

```
g chr3:187,597,641
```

Notice the \*s in the contig sequence, indicating a gap.

# Exercise 7: Browse the contig to reference alignments using IGV

Start IGV.

Select "View -> Preferences -> Alignments" and change "Visibility range threshold (kb)" to 1000. Select the "Genomes -> Load Genome From Server -> Human hg38". Then select "File -> Load from File" and "k48/HS0674-contigs.bam" Go to the region `chr3:186,500,000-188,000,000` by entering it into the box labeled "Go". IGV may take up to a minute to load.

Do the contigs overlap coding genes?

Add the dbSNP track. Select "File -> Load from Server", expand "Annotations" and select "All Snps 1.4.2". Should you have issues loading it, select "common Snps 1.4.2" instead.

Zoom in on a SNV. Is it in dbSNP?
