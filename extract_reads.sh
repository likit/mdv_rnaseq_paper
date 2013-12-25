#!/bin/sh

function usage {
    printf "Usage: %s <BAM file>\n" $0
    printf "A BAM file is required.\n"
} >&2

if [ $# -lt 1 ]
then
    usage
    exit 1
fi

if [ ! -e "$1".bai ]
then
    printf "BAM index file is not found, indexing it...\n" >&2
    samtools index "$1"
fi

outdir=$(dirname "$1")

if [ $# -eq 2 ]
then
    printf "A list of chromosome is provided in %s.\n" "$2" >&2
    for chr in $(cat "$2")
    do
        printf "\textracting reads for %s..." "$chr" >&2
        samtools view -o "$outdir"/"$chr".sam "$1" "$chr"
        printf "done.\n" >&2
    done
else
    samtools view -o "$outdir"/"reads".sam "$1"
fi
