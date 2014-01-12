'''Reads output from RSEM isoform expression and filter out
sequences with low isoform percentage.

'''

import sys
from collections import defaultdict

if len(sys.argv) < 3:
    print >> sys.stderr, \
        'Usage: filter-low-isopct.py cutoff bedfile isoforms-results'
    sys.exit(1)

cutoff = float(sys.argv[1])
bedfile = sys.argv[2]
db = defaultdict(bool)
for infile in sys.argv[3:]:
    fp = open(infile)
    _ = fp.readline()
    for line in fp:
        items = line.split('\t')
        transid = items[0]
        isopct = float(items[-1])
        if isopct > cutoff:
            db[transid] = True

for line in open(bedfile):
    seqid = line.split('\t')[3]
    if seqid in db:
        print line.strip()
