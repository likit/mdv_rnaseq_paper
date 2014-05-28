'''Parses Blast output in text format and extract

positions of SNPs between two sequences

'''

import sys

handle = open(sys.argv[1])

from Bio.Blast import NCBIStandalone
blast_parser = NCBIStandalone.BlastParser()
blast_record = blast_parser.parse(handle)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        i = 0
        while i in range(len(hsp.query)):
            try:
                q = hsp.query[i]
            except IndexError:
                print i, len(hsp.query)
                raise SystemExit

            if hsp.query[i] != hsp.sbjct[i]:
                j = i + 1
                while hsp.query[j] != hsp.sbjct[j]:
                    j += 1
                    if j > len(hsp.query) - 1:
                        break

                if j != (i + 1):  # run-on SNPs or indels
                    print '%d...%d\t%s\t%s' % (i + 1, j + 1,
                            hsp.sbjct[i:j], hsp.query[i:j])
                    i = j
                    continue
                else:  # a single SNP
                    # SNPs are in 1-based position
                    print '%d\t%s\t%s' % (i + 1,
                            hsp.sbjct[i], hsp.query[i])
                    i += 1
            else:  # no SNPs
                i += 1
