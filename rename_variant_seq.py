'''Renames variant sequences in FASTA format with IDs from GFF file'''

import sys
from csv import reader

ids = []

for rec in reader(open(sys.argv[1]), 'excel-tab'):
    i = rec[-1].split(';')[0].split('=')[1]
    ids.append(i)

j = 0
for line in open(sys.argv[2]):
    if line.startswith('>'):
        line = '>%s\n' % (ids[j])
        j += 1
    print line.strip()
