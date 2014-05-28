import sys
import csv
from collections import namedtuple

ProtDom = namedtuple('ProtDom', ['id', 'start', 'end',
                                'source', 'domainID', 'domain'])
colors = ['255,153,51', '0,204,0', '51,51,255',
        '255,51,153', '102,102,0', '102,51,0']
color = 0
source_colors = {}
domains = {}

with open(sys.argv[1]) as bedfile:
    reader = csv.reader(bedfile, dialect='excel-tab')
    for line in reader:
        # add 1 to start positions because BED is 0-based
        regid = '%s:%d-%s' % (line[0], int(line[1]) + 1, line[2])
        source = line[3]
        if source not in source_colors:
            try:
                source_colors[source] = colors[color]
                color += 1
            except IndexError:
                color = 0
                source_colors[source] = colors[color]
                color += 1

        domains[regid] = ProtDom(*line)

writer = csv.writer(sys.stdout, 'excel-tab')
with open(sys.argv[2]) as pslfile:
    reader = csv.reader(pslfile, dialect='excel-tab')
    for line in reader:
        id, start, end = line[9], int(line[15]), int(line[16])
        tbase_ins = int(line[7])

        block_counts = line[17].strip(',')
        chrom_size = int(line[14])

        # For dnax, each block size needs to be multiplied by three
        block_sizes = [int(i) * 3 for i in line[18].split(',')[:-1]]

        if line[8] == '+-':
            block_sizes.reverse()
            tstarts = [chrom_size - int(i) for i in line[20].split(',')[:-1]]
            tstarts = [i - start for i in tstarts]  # target start must be relative to start
            tstarts.reverse()  # target starts are backward, need to be reversed
            for i in range(len(block_sizes)):
                tstarts[i] = tstarts[i] - block_sizes[i]
        elif line[8] == '++':
            # target start must be relative to start
            tstarts = [int(i) - start for i in line[20].split(',')[:-1]]
            if tbase_ins < 0:  # deletion
                end += (tbase_ins * -1)
                tstarts[-1] += 1


        input = '\t'.join(line)
        error_msg = ' '.join(map(str, [block_sizes[-1], start, tstarts[-1], end]))
        assert int(block_sizes[-1]) + start + int(tstarts[-1]) == end, error_msg + '\n' + input  # sanity check

        tstarts = ','.join([str(i) for i in tstarts])
        block_sizes = ','.join([str(i) for i in block_sizes])

        strand = line[8][-1]
        target_name = line[13]

        if id in domains:
            dom = domains[id]
            chrom = dom.id.split(';')[0].split(':')[0]
            description = '%s:%s:%s' % (dom.source, dom.domainID, dom.domain)
            description = description.replace(' ', '_')
            if chrom != target_name:
                continue
            else:
                writer.writerow((chrom, start, end,
                            description, 1000, strand, start, end,
                            source_colors[dom.source],
                            block_counts, block_sizes, tstarts))
