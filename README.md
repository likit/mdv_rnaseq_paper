mdv_rnaseq_paper
================

Required software
----------------

+ Gimme v.0.98
+ CDHIT v.4.5.6
+ BLAT
+ Seqclean x86_64
+ Cufflinks v.2.1.1
+ Tophat2
+ Bowtie 1.0.0
+ RSEM 1.2.7
+ ESTScan 3.0.3
+ Seqtk
+ InterproScan
+ Velvet 1.2.03
+ Oases 0.2.06
+ BLAST+ 2.2.25
+ MISO 0.49
+ GATK 2.5.2
+ PicardTools 1.113
+ Samtools 0.1.19

#Protocol

###Setup paths

    export PROTOCOL=<path to the protocol>
    export GIMMEDIR=<path to Gimme>

At your working directory run all the following commands.

###Grab gene models and some prerequisite data

Please consult https://github.com/likit/RNASeq-methods-comparison
protocol on how to build Gimme models.

    wget http://athyra.ged.msu.edu/~preeyano/gene-network/pipeline/gimme/gimme.*
    wget http://athyra.ged.msu.edu/~preeyano/gene-network/pipeline/gal4selected*

###Prepare the reference genome

This step sorts the FASTA file and builds
a dictionary required for a SNP analysis using GATK.

    make -f $PROTOCOL/makefile protocol=$PROTOCOL init

###Prepare Reads

Get reads

    wget http://athyra.ged.msu.edu/~preeyano/mdv/paired-end/*
    wget http://athyra.ged.msu.edu/~preeyano/mdv/single-end/*

Quality trim

    make -f $PROTOCOL/misc.mk protocol=$PROTOCOL run-quality-trim-pe
    make -f $PROTOCOL/misc.mk protocol=$PROTOCOL run-quality-trim-se

###Run DE analysis

Build RSEM reference

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-prepare-ref

Run RSEM calculate expression

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-calc-expression

Run EBSeq

    make -f $PROTOCOL/makefile ebseq-line7 ebseq-line6

Convert RSEM DE output to FASTA

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-output-to-fasta

Get longest sequences

    make -f $PROTOCOL/makefile protocol=$PROTOCOL get-longest-sequences

Annotate sequences with chicken ENSEMBL genes

    make -f $PROTOCOL/makefile protocol=$PROTOCOL run-blast-gallus

###MISO analysis

Map reads to the genome

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL run-tophat-se

Merge BAM files from paired- and single-end reads and index them

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL merge-bam-files

Filter out low abundance isoforms

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL filter-low-isopct

Build SE, A3SS and A5SS models

    make -f $PROTOCOL/miso.mk gimmedir=$GIMMEDIR build-se-models build-a3ss-models build-a5ss-models

Run MISO

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL run-miso-se run-miso-a3ss run-miso-a5ss

Summarize MISO results

    make -f $PROTOCOL/miso.mk summarize-se summarize-a3ss summarize-a5ss

Compare MISO results

    make -f $PROTOCOL/miso.mk compare-miso-se compare-miso-a3ss compare-miso-a5ss

Filter MISO results

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL filter-miso-se filter-miso-a3ss filter-miso-a5ss

Annotate isoforms

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL miso-to-fa-se miso-to-fa-a3ss miso-to-fa-a5ss

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL blast-miso-se blast-miso-a3ss blast-miso-a5ss

Translate isoforms

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL translate-isoforms-miso-se translate-isoforms-miso-a3ss translate-isoforms-miso-a5ss

Run Interpro scan

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL interpro-isoforms-miso-se interpro-isoforms-miso-a3ss interpro-isoforms-miso-a5ss

BLAT domains

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL blat-domains-se blat-domains-a3ss blat-domains-a5ss

Annotate protein domains

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL annotate-domains-se annotate-domains-a3ss annotate-domains-a5ss

MISO to KEGG

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL miso-to-kegg-se miso-to-kegg-a3ss miso-to-kegg-a5ss

Find DEU snps

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL miso-snps-se miso-snps-a3ss miso-snps-a5ss
    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL projpath=$PWD find-deu-snps-se find-deu-snps-a3ss find-deu-snps-a5ss
