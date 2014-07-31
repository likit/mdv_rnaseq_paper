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

###Prepare RSEM reference

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-prepare-ref

###Run RSEM calculate expression

    make -f $PROTOCOL/makefile protocol=$PROTOCOL rsem-calculate-expression
