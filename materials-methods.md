#Materials and methods

##Sequences and quality trimming

mRNAs were extracted from spleens of control and infected
chickens lines 6 and 7 (4 d.p.i).  Sequence libraries were
prepared by standard Illumina unstranded single- and paired-end
protocols.  Library size of the paired-end datasets is
approximately 175 bp.  Read lengths are 75 bp in both single- and
paired-end libraries.  Reads were quality trimmed by Condetri 2.1
with quality score cutoff of 30.  The first 10 bases were removed
due to non-uniform distribution of nucleotides.

##Gene model construction

Due to lack of complete gene models for chickens, we employed two
methods to construct gene models from RNA-Seq reads.  First,
short reads were assembled using Velvet/1.2.03 and Oases/0.2.06
to obtain long transcripts.  Assembly was done with hash lengths
range from 21 to 31 for both local and global assembly (described
in Gimme paper).  Poly-A tails were trimmed and low complexity
transcripts were removed by Seqclean.  All transcripts were then
aligned to the chicken reference genome (galGal4, with unplaced
and random chromosomes removed) with BLAT.  Second, reads were
aligned to the reference genome using Tophat/2.0.9 and ENSEMBL
gene model release 73 was used to guide reference-based assembly
by Cufflinks2.  Alignments from BLAT and models from Cufflinks
were then combined to construct gene models by Gimme (manuscript
in preparation).

##Differential gene expression analysis and gene ontology

To identify DE genes, reads were mapped to transcripts by RSEM
v.1.2.7, which is also used to estimate gene expression and
identify DE genes.  Data from single- and paired-end datasets
from the same line were treated as biological replicates.  To
identify enriched pathways and ontology terms, a list of DE genes
was analysed by GOSeq v.1.10.0 based on chicken KEGG annotations.
P-values were corrected by Benjamini-Hochberg multiple testing
correction.  Genes, pathways and GO terms with corrected
P-value<0.1 were considered significant.  Pathview was used to
create a KEGG pathway diagram with colors representing relative
level of gene expressions.

##Differential exon usage analysis

Gene models were converted to alternative splicing models using a
Python script.  In order to increase sensitivity, read counts
from single- and paired-end samples were combined and treated as
single-end reads for splicing event analysis with MISO/0.4.9.
Splicing events with Bayes factor $>10$ and delta Psi>0.20 were
considered significant.  Read coverages and Psi distributions
were plotted using Sashimi plot.

##Variant calling and __in silico__ splicing analysis

Variants were called using mpileup command from
SAMTools/0.1.18 and
BCFTools.  Only variants with quality score >= 20 were used for
mutation analyses.  Exon enhancers and suppressors were predicted
using the Human Splicing Finder web portal.  Human default
parameter settings were used in all analyses. 

##Protein domains search

Transcripts were translated to protein sequences using ESTScan
3.0.3 with chicken HMM matrices built from chicken cDNAs and
refSeq sequences. Protein sequences were searched against
InterPro database using InterProScan/5.44.0.
