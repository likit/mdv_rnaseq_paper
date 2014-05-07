import sys
from Bio import SeqIO

genes = set()
for line in open(sys.argv[2]):
    genes.add(line.strip())

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    desc = rec.description.split()
    gene, transcript = desc[3].split(':')[1], desc[4].split(':')[1]
    if gene in genes:
        rec.id = '%s:%s' % (gene, transcript)
        rec.description = ''
        SeqIO.write(rec, sys.stdout, 'fasta')
