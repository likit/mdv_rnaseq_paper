import sys

for line in sys.stdin:
    items = line.strip().split('\t')
    chr = items[0].lstrip('chr')
    pos = items[1]
    ref = items[3]
    alt = items[4]
    line7 = items[7] if items[7] != ref else "."
    strand = items[-3]
    eventid = items[-1].split(';')[0].split('=')[1]
    print ','.join([eventid, chr, pos, ref, alt, line7, strand])
