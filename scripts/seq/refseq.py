#!/usr/bin/env python2.7

# from Bio import SeqIO
from collections import defaultdict

fRefseq = open('refseq.fasta')
fRef2symbol = open('ref2symbol.txt')
# for seq in SeqIO.parse(f, "fasta"):
#     print(record.id, len(record))

# f.close()

dRef2symbol = defaultdict(str)

for line in fRef2symbol.readlines():
	line = line.rstrip('\n')
	line = line.split('\t')
	dRef2symbol[line[0]] = line[1]

for line in fRefseq.readlines():
	line = line.rstrip('\n')
	if line.startswith('>'):
		refId = line.split('|')[3]
		refId = refId.split('.')[0]
		print '>' + dRef2symbol[refId], refId
	else:
		print line

fRefseq.close()
fRef2symbol.close()