# Iso-seq data pipeline

## PacBio concepts
PacBio SMRT sequencing operates within a silicon chip (**SMRTcell**) fabricated to contain a large number of microscopic holes (**ZMWs**, or zero-mode waveguides), each assigned a hole number.

Within a ZMW, PacBio SMRT sequencing is performed on a circularized molecule called a **SMRTbell**. The SMRTbell, depicted below, consists of:

* the customer’s double-stranded DNA insert (I)
* (optional) double-stranded DNA barcodes (sequences BL, BR) used for multiplexing DNA samples. While the barcodes are optional, they must be present at both ends if present at all.
* SMRTbell adapters (sequences AL, AR), each consisting of a double stranded stem and a single-stranded hairpin loop. 
  
![image](https://user-images.githubusercontent.com/31697487/135611262-770421a4-0538-4738-865b-750f25bae351.png)
  
Specifically the raw data or each SMRT-cell will be in files named  ```.subreads.bam```,  ```.subreads.xml```, and  ```.subreads.pbi```.

The ```.bam data``` can be converted to fastq or fasta files with bamtools.


### What is Circular Consensus Sequence or ccs?
ccs combines multiple subreads of the same SMRTbell molecule using a statistical model to produce one highly accurate consensus sequence, also called a HiFi read.
Single Molecule, Real-Time (SMRT) Sequencing technology has evolved to a different type of long read, known as highly accurate long reads, or **HiFi reads**. PacBio is the only sequencing technology to offer HiFi reads that provide accuracy of >99.9%, on par with short reads and Sanger sequencing.
 
 
![alt text](https://ccs.how/img/generate-hifi.png)


## Step 0: Input
For each SMRT cell a ```duck.subreads.bam``` is needed for processing.

## Step 1: Circular Consensus Sequence (ccs) calling 

Each sequencing run is processed by ccs to generate one representative circular consensus sequence (CCS) for each ZMW.
* **Input**: Subreads from a single movie in PacBio BAM format (```.subreads.bam```).
* **Output**: Consensus reads in a format inferred from the file extension: unaligned BAM (```.bam```); bgzipped FASTQ (```.fastq.gz```); 
or SMRT Link XML (```.consensusreadset.xml```) which also generates a corresponding unaligned BAM file.

more info at [how does ccs work](https://ccs.how/how-does-ccs-work.html)
```bash
ccs duck.subreads.bam duck.ccs.bam --min-rq 0.9
```
### Parallelize by chunk

For this, the ```.subreads.bam``` file must accompanied by a ```.pbi``` file. To generate the index ```subreads.bam.pbi```, use ```pbindex```, which can be installed with ```conda install pbbam```.
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
Output files will be called according to their primer pair. Example for single sample libraries: ```duck.fl.NEB_5p--NEB_Clontech_3p.bam```

## Step 3a: Refine
Data now contains full-length (FL) reads, but still needs to be refined by:
* Trimming of poly(A) tails.
* Rapid concatemer identification and removal.

So, input and output for this step should be:
* **Input**: The input file for refine is one demultiplexed CCS file with full-length reads and the primer fasta file.
* **Output**: The output files of refine contain *Full Length Non-Chimeric (FLNC)* reads.
```bash
isoseq refine duck.NEB_5p--NEB_Clontech_3p.fl.bam primers.fasta duck.flnc.bam
```

If your sample has poly(A) tails, use ```--require-polya```. This filters for FL reads that have a poly(A) tail with at least 20 base pairs (```--min-polya-length```) and removes identified tail:
```bash
isoseq refine duck.NEB_5p--NEB_Clontech_3p.fl.bamduck.flnc.bam --require-polya
```

## Step 3b: Merge SMRT Cells (optional)
If you used more than one SMRT cells, list all of your ```<duck>.flnc.bam``` in one ```flnc.fofn```, a file of filenames:
```bash
ls duck*.flnc.bam duck*.flnc.bam duck*.flnc.bam > flnc.fofn
 ```

## Step 4: Clustering (optional)
This step takes the FLNC reads and clusters the transcript sequences by similarity. It then makes a multiple alignment of each cluster and performs error correction using this alignment. This step is optional and if you are interested in allele specific expression/transcript phasing, you should not run this step as it can removed allele specific sequence variation. Give this step as many cores as possible

* **Input** The input file for cluster is one FLNC file:
```<duck>.flnc.bam``` or ```flnc.fofn```

* **Output** The following output files of cluster contain polished isoforms:

   *  ```duck.bam ```
   *  ```duck.hq.fasta.gz ```  with predicted accuracy ≥ 0.99
   *  ```duck.lq.fasta.gz ```  with predicted accuracy < 0.99
   *  ```duck.bam.pbi ```
   *  ```duck.transcriptset.xml ```
 
```bash
isoseq cluster flnc.fofn clustered.bam --verbose --use-qvs
  ```
   
more info of [Isoseq3 pipeline](https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md)

## Step 4: Mapping reads against genome with Minimap2
This step can take inputs from either the FLNC reads or the Polished reads. Minimap2 maps the transcript sequences onto a genome assembly.
```bash
minimap2 -ax splice -uf -C5 genome_duck_ref.fa duck_iso-seq.fq > duck_aln.sam
```
Option ```-C5``` reduces the penalty on non-canonical splicing sites. It helps to align such sites correctly for data with low error rate such as Iso-seq reads and traditional cDNAs. On this example, minimap2 makes one junction error. Applying ```--splice-flank=no``` fixes this alignment error.
splice:hq?
  
more info in [Minimap2](https://github.com/lh3/minimap2)

Note that the command line above is optimized for the final Iso-seq reads. PacBio's Iso-seq pipeline produces intermediate sequences at varying quality. For example, some intermediate reads are not stranded. For these reads, option -uf will lead to more errors. Please revise the minimap2 command line accordingly.

## Step 5: Collapse 
TAMA Collapse is a tool that allows you to collapse redundant transcript models in your Iso-Seq data.
This step takes a bam/sam file from the transcript mapping and collapses redundant transcripts based on genomic location.

Also provides the following information: 
* Source information for each predicted feature
* Variation calling
* Genomic poly-A detection 
* Strand ambiguity

```bash
python $tama_dir/tama_collapse.py -s $tamaout_path/$name.aln.sort.sam -f $genome_ref -p $name_out -x capped -a 100 -z 100 -lde 5 -sj sj_priority -sjt 20
```
**x:** If you performed 5' cap selection during your Iso-Seq library preparation you should use the "-x capped" option. This will not allow transcripts to be collapses if they have a different number of exons.  
**a:** The 5 prime threshold is the amount of tolerance at the 5' end of the transcript for grouping reads to be collapsed.  
**z:** The 3 prime threshold is the amount of tolerance for the 3' end of the transcript for grouping reads to be collapsed.  
**sj:** I call this option splice junction priority. Basically if you turn this on, TAMA Collapse will search for mismatches in the regions on each side of each splice junction. The length of the regions are defined by the SJT flag. TAMA Collapse will then rank the splice junction evidence from each transcript/read as either 0, 1, 2, or 3 with 0 being highest priroity.   
**lde:** This is the Threshold for the number of errors that is allowed to occur on each side of the splice junction within the range specified by the -sjt argument.This feature is used to remove mapped reads where the splice junction mapping is likely to be erroneous due to a high local density of errors near the splice junction.  
**sjt:** This is the length threshold for the regions to observe for splice junction priority. The default is 10bp which means TAMA collapse will look at the 10bp region upstream and downstream of the splice junction.  

Command used in human transcriptome: TAMA collapse (−d merge_dup -x no_cap -a 100 -z 100 -sj sj_priority -lde 5 -sjt 20 -log log_off)
more info of [TAMA Collapse](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse)

## Step 5: Merge
Merging is the process of combining multiple transcriptomes. For instance, if you have Iso-Seq data for different tissue types you might want to process them separately and then combine them at the end to use as a transcriptome for downstream analysis. However, the act of merging transcriptomes is non-trivial with respect to choosing what criteria to use to merge transcript models. You probably would also like to keep a record of which models from your merged transcriptome came from which source. 
```bash
python tama_merge.py -f filelist.txt -p merged_annos
```
NOTE: If you do not see "TAMA Merge has completed successfully!" as the last line of the terminal output, then TAMA Merge has not completed and there are likely issues with the input files.

more info of [TAMA Merge](https://github.com/GenomeRIK/tama/wiki/Tama-Merge)

## Step 6: TAMA GO
### 1. Converting the bed files into fasta files
Use Bedtools to convert your bed file into a fasta file. (Bedtools 2.26.0 or newer)

```bash
${bedtools} getfasta -name -split -s -fi ${genome_ref} -bed ${bed} -fo ${outfile}
```
fasta - this is the fasta file for the genome assembly  
bed - this is the bed12 file for the annotation  
outfile - this is the name of the output fasta file  

The output fasta headers should look like this:

```bash
  >G1;G1.1::1:29-2665(+)
  >G1;G1.2::1:219-3261(+)
  >G1;G1.3::1:1713-3246(+)
```
### 2. Getting open read frames (ORF) from the transcript sequences
Use tama_orf_seeker.py to get ORF amino acid sequences.

```bash
python ${tama_dir}/tama_go/orf_nmd_predictions/tama_orf_seeker.py -f ${mergefasta} -o ${outfile}_orf
```

fasta - this is the fasta file for your transcript sequences that you produced in the previous step  
outfile - the name of the output fasta file which contains the ORF amino acid sequences  

The output fasta file should have headers like this:
```bash
  >G1;G1.1::1:29-2665(+):F1:0:719:1:239:239:I
  >G1;G1.1::1:29-2665(+):F1:132:719:45:239:195:M
  >G1;G1.1::1:29-2665(+):F3:2:64:1:20:20:F
```

### 3. Blasting the amino acid sequences against the Uniprot/Uniref protein database
Use blastp and a local Uniref database to find protein hits. (BLAST 2.2.31+ or newer)  
uniref - Uniref/Uniprot fasta file  
orf - ORF fasta file from previous step  
```bash
${blastp} -evalue 1e-10 -num_threads 16 -db ${uniref} -query ${orf} > ${name_out}.blastp 
```

The individual hit entries in the output file should look like this:
```bash
  >E1C721 Uncharacterized protein OS=Gallus gallus GN=LOC425783 PE=4 SV=2
  Length=265
  
   Score = 382 bits (982),  Expect = 4e-136, Method: Compositional matrix adjust.
   Identities = 203/235 (86%), Positives = 209/235 (89%), Gaps = 12/235 (5%)
    
  Query  17   KSVLAKKSAPPAPLCPQPGPSLPLSPHTMGAVPHLHGAKGERERRSPSPPREATAREGDE  76
              + VLAKKSAPPAPLCPQPGPSLPLSPHTMGAVPHLHGAKGERERRSPSPPREATAREGDE
  Sbjct  31   REVLAKKSAPPAPLCPQPGPSLPLSPHTMGAVPHLHGAKGERERRSPSPPREATAREGDE  90
```
### 4. Parsing the Blastp output file for top hits
Use tama_orf_blastp_parser.py to parse the results from Blastp.

```bash
python ${tama_dir}/tama_go/orf_nmd_predictions/tama_orf_blastp_parser.py -b ${name_out}.blastp -o ${outfile}_blastparser
```

If you did a BlastP on an Ensembl CDS peptide sequence database, you can add "-f ensembl" to the command line to allow for parsing of the Ensembl ID's.

The output should look like this:  
```bash
  G1;G1.1 F1      0       719     1       239     R4GL70  90_match        99      32
  G1;G1.2 F1      369     596     124     198     R4GIW7  50_match        77      66
  G1;G1.3 F1      0       386     1       128     H9KZN5  90_match        91      65
  ```
### 5. Create new bed file with CDS regions
Use tama_cds_regions_bed_add.py to create the annotation bed file with CDS regions added.  
orf - The output file from the previous step  
bed - The bed annotation file that you started with in this pipeline  
fasta - The transcript input fasta file from the first step (the input fasta)  
outfile - The output bed file with CDS regions  

```bash
python ${tama_dir}/tama_go/orf_nmd_predictions/tama_cds_regions_bed_add.py -p ${outfile}_blastparser -a ./${name_out}/${bed} -f ${mergefasta} -o ${name_out}.cds
```

### 6. Explanation of final output
The output of step 5 is a bed file with CDS regions defined by column 7 and 8. Column 7 is the start of the CDS region and column 8 is the end. The 4th column contains the gene ID and transcript ID along with additional information.

This is the format of the 4th column:  

**gene_id;trans_id;gene_name;trans_length_flag;blastp_match_flag;nmd_flag;frame**

This is an example of a line:
```bash
  1       219     3261    G1;G1.2;R4GIW7;full_length;50_match;prot_ok;F1  40      +       2154    2662    200,0,255       5       98,93,181,107,714       0,1457,1757,2132,2328
```
