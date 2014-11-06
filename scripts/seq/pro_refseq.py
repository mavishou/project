#!/usr/bin/env python2.7

# from Bio import SeqIO
from collections import defaultdict

fRefseq = open('pro_refseq.fa')
fSymbol2ref = open('pro_symbol2ref.txt')
# for seq in SeqIO.parse(f, "fasta"):
#     print(record.id, len(record))

# f.close()

dRef2symbol = defaultdict(str)

for line in fSymbol2ref.readlines():
	line = line.rstrip('\n')
	line = line.split('\t')
	dRef2symbol[line[1]] = line[0]

for line in fRefseq.readlines():
	line = line.rstrip('\n')
	if line.startswith('>'):
		refId = line.split('|')[3]
		refId = refId.split('.')[0]
		print '>' + dRef2symbol[refId], refId
	else:
		print line

fRefseq.close()
fSymbol2ref.close()