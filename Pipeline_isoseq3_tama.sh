

#### Step 1 - Circular Consensus Sequence calling
#Each sequencing run is processed by ccs to generate one representative circular consensus sequence (CCS) for each ZMW.

ccs movieX.subreads.bam movieX.ccs.bam --min-rq 0.9

#### Step 2 - Primer removal and demultiplexing
#Removal of primers and identification of barcodes is performed using lima that offers a specialized --isoseq mode. 
#If there are more than two sequences in your primer.fasta file or better said more than one pair of 5' and 3' primers.
#Use lima with --peek-guess to remove spurious false positive signal.

lima movieX.ccs.bam barcoded_primers.fasta movieX.fl.bam --isoseq --peek-guess

Lima will remove unwanted combinations and orient sequences to 5' → 3' orientation.

Output files will be called according to their primer pair. Example for single sample libraries: movieX.fl.NEB_5p--NEB_Clontech_3p.bam

#If your library contains multiple samples, execute the following workflow for each primer pair:

movieX.fl.primer_5p--brain_3p.bam
movieX.fl.primer_5p--liver_3p.bam

#### Step 3 - Refine
# data now contains full-length reads, but still needs to be refined by:
# a)Trimming of poly(A) tails
# b)Rapid concatemer identification and removal
# The input file for refine is one demultiplexed CCS file with full-length reads and the primer fasta file:

<movie.primer--pair>.fl.bam or <movie.primer--pair>.fl.consensusreadset.xml
primers.fasta
Output The following output files of refine contain full-length non-concatemer reads:

<movie>.flnc.bam
<movie>.flnc.transcriptset.xml
Actual command to refine:

$ isoseq refine movieX.NEB_5p--NEB_Clontech_3p.fl.bam primers.fasta movieX.flnc.bam
If your sample has poly(A) tails, use --require-polya. This filters for FL reads that have a poly(A) tail with at least 20 base pairs (--min-polya-length) and removes identified tail:

$ isoseq refine movieX.NEB_5p--NEB_Clontech_3p.fl.bam movieX.flnc.bam --require-polya
Step 3b - Merge SMRT Cells
If you used more than one SMRT cells, list all of your <movie>.flnc.bam in one flnc.fofn, a file of filenames:

$ ls movie*.flnc.bam movie*.flnc.bam movie*.flnc.bam > flnc.fofn
