import sys
import csv

try:
    inputfile = sys.argv[1]
except IndexError:
    print >> sys.stderr, 'No input file found, reading data from stdin.'
    reader = csv.reader(sys.stdin, dialect="excel-tab")
else:
    reader = csv.reader(open(inputfile), dialect="excel-tab")

writer = csv.writer(sys.stdout, dialect="excel-tab")
for line in reader:
    writer.writerow([line[0], line[6], line[7], line[3], line[4], line[5]])
