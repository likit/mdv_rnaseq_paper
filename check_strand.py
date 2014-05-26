import sys
import csv
from Bio import SeqIO

reader = csv.reader(open(sys.argv[1]), 'excel-tab')
prefix = sys.argv[2]

for rec in reader:
    strand = rec[6]
    gene_id = rec[-1].split(';')[0].split('=')[1]
    # print gene_id, strand
    seq = SeqIO.read('%s.%s.fa' % (prefix, gene_id), 'fasta')

    if strand == '-':
        seq = seq.reverse_complement()

    seq.id = gene_id
    seq.description = strand
    SeqIO.write(seq, '%s.%s.chk.fa' %
                (prefix, gene_id), 'fasta')
