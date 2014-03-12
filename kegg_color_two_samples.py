import sys

sample1 = sys.argv[1]
sample2 = sys.argv[2]

data1 = set()
data2 = set()

for line in open(sample1):
    geneid = line.split('\t')[1]
    data1.add(geneid)

for line in open(sample2):
    geneid = line.split('\t')[1]
    data2.add(geneid)

for geneid in data2.intersection(data1):
    print '%s\tgreen' % (geneid)

for geneid in data1.difference(data2):
    print '%s\tyellow' % (geneid)

for geneid in data2.difference(data1):
    print '%s\tblue' % (geneid)

