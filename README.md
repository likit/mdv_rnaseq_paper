mdv_rnaseq_paper
================

Required software
----------------

+ Gimme v.0.98
+ CDHIT v.4.5.6
+ BLAT
+ Seqclean x86_64
+ Cufflinks v.2.1.1

#Protocol

###Setup paths

    export PROTOCOL=<path to the protocol>
    export GIMMEDIR=<path to Gimme>

At your working directory run all the following commands.

###Prepare Reads

    make -f $PROTOCOL/misc.mk protocol=$PROTOCOL run-quality-trim-pe
    make -f $PROTOCOL/misc.mk protocol=$PROTOCOL run-quality-trim-se

###Cufflinks

Map reads to the genome using Tophat

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-tophat-pe
    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-tophat-se

Run RSEM prepare reference

    make -f $PROTOCOL/cufflinks.mk prepare-reference-cufflinks

Then run Cufflinks

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-cufflinks

Merge gene models from all samples

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-cuffmerge-ref

###Local assembly

Extract reads from each chromosome

    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL extract-reads

Merge reads from each chromosome to local assembly directory

    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL merge-bams

Run Velveth, Velvetg and Oases

    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL run-velveth-local
    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL run-velvetg-local
    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL run-oases-local

Merge transcripts

    make -f $PROTOCOL/local_assembly.mk gimmedir=$GIMMEDIR combine-transcripts

###Merged models

Combine all transcripts

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL gimmedir=$GIMMEDIR combine-transcripts

Clean transcripts

    make -f $PROTOCOL/gimme.mk clean-transcripts

Remove redundant sequences

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL remove-redundant-seq

Align transcripts

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL align-transcripts

Sort and select best alignments

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL sort-alignments

Build merged models

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL gimmedir=$GIMMEDIR build-merged-gene-models

Note, gimme.bed may contain warning messages from pygr at the beginning of the file.
The messages have to be removed before running the next step.

Get transcripts from gene models

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL gimmedir=$GIMMEDIR merged-models-to-transcripts

###Run DE analysis

Build RSEM reference

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-prepare-ref

Run RSEM calculate expression

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-calculate-expression
