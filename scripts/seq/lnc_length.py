#!/usr/bin/env python2.7

from Bio import SeqIO

fLnc = open('lncRNA_final.fa')
fRPI = open('RPI_final.txt')
dLncLen = {}

for record in SeqIO.parse(fLnc, "fasta"):
	dLncLen[record.id] = len(record)

for line in fRPI.readlines():
	line = line.rstrip('\n').split('\t')
	out = '\t'.join([line[0], line[1], str(dLncLen[line[0]])])
	print out

fLnc.close()
fRPI.close()