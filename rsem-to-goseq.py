'''Create a genes list for GOSeq analysis from RSEM output.'''

import sys

if len(sys.argv) < 3:
    print >> sys.stderr, 'Usage: rsem-to-goseq.py bedfile rsem-output'
    sys.exit(1)

bedfile = sys.argv[1]
rsemout = open(sys.argv[2])
degenes = set()
allgenes = set()

_ = rsemout.readline()  # skip the header

for line in rsemout:
    geneid = line.split('\t')[0].replace('"', '')
    degenes.add(geneid)

for line in open(bedfile):
    geneid  = line.split('\t')[3].split('.')[0]
    allgenes.add(geneid)

for geneid in allgenes:
    if geneid in degenes:
        print '%s\t1' % geneid
    else:
        print '%s\t0' % geneid
