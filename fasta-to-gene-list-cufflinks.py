import sys

for line in sys.stdin:
    if line.startswith('>'):
        tranid = line.lstrip('>').strip()
        geneid = '.'.join(tranid.split('.')[:2])
        print '%s\t%s' % (geneid, tranid)
