#!/usr/bin/env python
'''Description'''
import sys
from Bio import SeqIO

def read_list(infile):
    seqlist = []

    for line in open(infile):
        seqlist.append(line.strip())

    return seqlist


def sort(infile, seqlist, outfile):

    seqdict = {}
    for rec in SeqIO.parse(infile, 'fasta'):
        seqdict[rec.id] = rec

    sortedseqs = []
    for seq in seqlist:
        sortedseqs.append(seqdict[seq])

    SeqIO.write(sortedseqs, outfile, 'fasta')


def main():
    '''Main function'''

    if len(sys.argv) < 4:
        print >> sys.stderr, \
            'Usage: sort-fasta.py <fasta> <list> <output filename>'
        raise SystemExit

    seqlist = read_list(sys.argv[2])
    outfile = sys.argv[3]
    sort(sys.argv[1], seqlist, outfile)


if __name__=='__main__':
    main()

