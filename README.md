mdv_rnaseq_paper
================

Required software
----------------

+ Gimme v.0.98
+ CDHIT v.4.5.6
+ BLAT
+ Seqclean x86_64
+ Cufflinks v.2.1.1
+ Tophat 2.0.9
+ Bowtie 2.2.1
+ RSEM 1.2.7
+ Bowtie 1.0.0
+ ESTScan 3.0.3
+ Seqtk
+ InterproScan 5.44.0
+ Velvet 1.2.03
+ Oases 0.2.06
+ BLAST+ 2.2.25
+ MISO 0.4.9
+ GATK 2.5.2
+ PicardTools 1.113
+ Samtools 0.1.19

#Protocol

This protocol has been tested on MSU HPC computer cluster only.
More info about the system can be found at https://icer.msu.edu/hpcc.

###Download scripts and Gimme

    git clone https://github.com/likit/mdv_rnaseq_paper.git mdv-protocol
    wget https://github.com/likit/gimme/archive/v.0.98.tar.gz
    tar xvfz v.0.98.tar.gz
    mv v.0.98 gimme

###Setup paths

    # please change these paths accordingly
    export PROTOCOL=mdv-protocol
    export GIMMEDIR=gimme

Create a working directory and run all the following commands inside the
directory.

    # please change the path accordingly
    mkdir mdv-project
    cd mdv-project

###Acquire gene models and some prerequisite data

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

Create Gallus gallus BLAST database:

    make -f $PROTOCOL/makefile create-gallus-blastdb

Annotate sequences with chicken ENSEMBL genes

    make -f $PROTOCOL/makefile protocol=$PROTOCOL projdir=$PWD run-blast-gallus
    make -f $PROTOCOL/makefile protocol=$PROTOCOL projdir=$PWD get-tophits

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

Find DEU snps

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL miso-snps-se miso-snps-a3ss miso-snps-a5ss
    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL projpath=$PWD find-deu-snps-se find-deu-snps-a3ss find-deu-snps-a5ss

Annotate-miso:

    make -f $PROTOCOL/miso.mk protocol=$PROTOCOL annotate-miso-se annotate-miso-a3ss annotate-miso-a5ss

Plot ITGB2 and PFN2 genes.

    plot.py --plot-event=chr7-24308.ev1 miso/indexes/SE \
    ~/mdv-protocol/sashimi_plot_settings_SE.txt \
    --plot-title=ITGB2 --output-dir=miso/results/plots

    plot.py --plot-event=chr9-25729.ev1 miso/indexes/A3SS \
    ~/mdv-protocol/sashimi_plot_settings_A3SS.txt \
    --plot-title=PFN2 --output-dir=miso/results/plots

###Run the actual analysis

Copy data

    # on your local computer
    mkdir -p mdvproj/results
    cd mdvproj; make -f $PROTOCOL/copy.mk all

Get gene info from chicken and human annotation

    cd mdvproj; make -f $PROTOCOL/analysis.mk annotate
    cd mdvproj; make -f $PROTOCOL/analysis.mk run-kegg

Now you should launch the IPython notebook server and interactively run the all
cells in the notebook.
