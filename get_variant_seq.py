'''The script reads features (gene, exon, etc) in GFF file

and converts them to an intervals in the format of
chr:start-end.

'''

import sys
import subprocess

from csv import reader, writer

if __name__=='__main__':

    gff_file = sys.argv[1]
    variant = sys.argv[2]
    outprefix = sys.argv[3]
    try:
        r = reader(open(gff_file), 'excel-tab')
    except:
        print >> sys.stderr, 'Cannot open the input file.'
        raise SystemExit

    for rec in r:
        id = rec[-1].split(';')[0].split('=')[1]
        interval =  "%s:%s-%s" % (rec[0], rec[3], rec[4])
        command = "GATK -T FastaAlternateReferenceMaker " \
                    "-R /mnt/ls12/preeyanon/mdv-pipeline/gal4selected.sorted.fa " \
                    "-o %s -V:testsnp,VCF %s -L %s" \
                    % (outprefix + '.' + id + '.fa', variant, interval)

        subprocess.call(command, shell=True)
