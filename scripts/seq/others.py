#!/usr/bin/euv python2.7

from Bio import SeqIO
from math import ceil

fOthers = open('others.fa')

def outputSeq(s, n = 70):
	nLine = int(ceil(float(len(s)) / n))
	for i in range(nLine):
		print s[(i * n) : (i + 1) * n]
	if (i + 1) * n < len(s):
		print s[(i + 1) * n : len(s)]

for record in SeqIO.parse(fOthers, "fasta"):
	print '>' + record.description
	seq = str(record.seq)
	seq = seq.replace(' ', '')
	seq = seq.replace('-', '')
	seq = seq.replace('U', 'T')
	outputSeq(seq)
	print
