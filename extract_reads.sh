#!/bin/sh

if [ $# -lt 1 ]
then
    printf "Usage: extract_reads.sh <BAM file>\n"
    printf "A BAM file is not required.\n"
    exit 1
fi

if [ ! -e "$1".bai ]
then
    printf "BAM index file is not found, indexing it...\n"
    samtools index "$1"
fi

outdir=$(dirname "$1")

if [ $# -eq 2 ]
then
    printf "A list of chromosome is provided in %s.\n" "$2"
    for chr in $(cat "$2")
    do
        printf "\textracting reads for %s..." "$chr"
        samtools view -o "$outdir"/"$chr".sam "$1" "$chr"
        printf "done.\n"
    done
else
    samtools view -o "$outdir"/"reads".sam "$1"
fi
