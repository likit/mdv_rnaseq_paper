import sys
from Bio import SeqIO

seqs = []
for n, seq in enumerate(SeqIO.parse(sys.argv[1], 'fasta')):
    if n % 10000 == 0:
        SeqIO.write(seqs, 'subsets_%d.fa' % n, 'fasta')
        seqs = []
    seqs.append(seq)

if seqs:
    SeqIO.write(seqs, 'subsets_%d.fa' % (len(seqs) + n), 'fasta')
