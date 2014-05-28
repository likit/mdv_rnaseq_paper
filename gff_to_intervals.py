'''The script reads features (gene, exon, etc) in GFF file

and converts them to an intervals in the format of
chr:start-end.

'''

import sys
from csv import reader, writer

try:
    r = reader(open(sys.argv[1]), 'excel-tab')
except:
    print >> sys.stderr, 'Cannot open the input file.'
    raise SystemExit

for rec in r:
    print '%s:%s-%s' % (rec[0], rec[3], rec[4])
