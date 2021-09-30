# Isoseq3 and TAMA pipeline with Iso-seq data

## Pacbio essencials

### What is Circular Consensus Sequence or ccs?
ccs combines multiple subreads of the same SMRTbell molecule using a statistical model to produce one highly accurate consensus sequence, also called a HiFi read.
Single Molecule, Real-Time (SMRT) Sequencing technology has evolved to a different type of long read, known as highly accurate long reads, or HiFi reads. PacBio is the only sequencing technology to offer HiFi reads that provide accuracy of >99.9%, on par with short reads and Sanger sequencing.
 
 
![alt text](https://ccs.how/img/generate-hifi.png)



## Step 1: Circular Consensus Sequence (ccs) calling

Each sequencing run is processed by ccs to generate one representative circular consensus sequence (CCS) for each ZMW.
Input: Subreads from a single movie in PacBio BAM format (```.subreads.bam```).
Output: Consensus reads in a format inferred from the file extension: unaligned BAM (```.bam```); bgzipped FASTQ (```.fastq.gz```); 
or SMRT Link XML (```.consensusreadset.xml```) which also generates a corresponding unaligned BAM file.

```bash
ccs duck.subreads.bam duck.ccs.bam --min-rq 0.9
```
### Parallelize by chunk

For this, the .subreads.bam file must accompanied by a ```.pbi``` file. To generate the index ```subreads.bam.pbi```, use ```pbindex```, which can be installed with ```conda install pbbam```.
```bash
pbindex duck.subreads.bam
```
An example workflow, all ccs invocations can run simultaneously:
```bash
ccs duck.subreads.bam duck.ccs.1.bam --chunk 1/10 -j <THREADS>
ccs duck.subreads.bam duck.ccs.2.bam --chunk 2/10 -j <THREADS>
...
ccs duck.subreads.bam duck.ccs.10.bam --chunk 10/10 -j <THREADS>
```
Merge chunks with pbmerge and index with pbindex
```bash
pbmerge -o duck.ccs.bam duck.ccs.*.bam
pbindex duck.ccs.bam
 ```
or use samtools
```bash
samtools merge -@8 duck.ccs.bam duck.ccs.*.bam
```
## Step 2: Primer removal and demultiplexing
Removal of primers and identification of barcodes is performed using lima that offers a specialized ```--isoseq``` mode. 
If there are more than two sequences in your primer.fasta file or better said more than one pair of 5' and 3' primers.
Use lima with ```--peek-guess``` to remove spurious false positive signal. Lima will remove unwanted combinations and orient sequences to 5' → 3' orientation.
```bash
lima duck.ccs.bam barcoded_primers.fasta duck.fl.bam --isoseq --peek-guess
```
Output files will be called according to their primer pair. Example for single sample libraries: ```movieX.fl.NEB_5p--NEB_Clontech_3p.bam```

## Step 3: Refine
data now contains full-length reads, but still needs to be refined by:

a)Trimming of poly(A) tails
b)Rapid concatemer identification and removal

The input file for refine is one demultiplexed CCS file with full-length reads and the primer fasta file:
The output files of refine contain full-length non-concatemer reads.
```bash
isoseq refine duck.NEB_5p--NEB_Clontech_3p.fl.bam primers.fasta duck.flnc.bam
```

If your sample has poly(A) tails, use --require-polya. This filters for FL reads that have a poly(A) tail with at least 20 base pairs (```--min-polya-length```) and removes identified tail:
```bash
isoseq refine duck.NEB_5p--NEB_Clontech_3p.fl.bamduck.flnc.bam --require-polya
```
