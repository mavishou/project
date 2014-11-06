#!/usr/bin/env python2.7
'''This script is used to tidy all protein sequences downloaded from ensembl.
Current filters:
1. sequence no *
2. sequence >50pp'''

from Bio import SeqIO
from math import ceil

cutoff = 50

def outputSeq(record, wOut, n = 70):
	'''a function for output sequence in correct format
	record: a seqIO object
	wOut: the output file
	'''
	s = str(record.seq)
	id = record.id
	wOut.write('>' + id + '\n')
	nLine = int(ceil(float(len(s)) / n))
	for i in range(nLine):
		# print s[(i * n) : (i + 1) * n]
		wOut.write(s[(i * n) : (i + 1) * n] + '\n')
	if (i + 1) * n < len(s):
		# print s[(i + 1) * n : len(s)]
		wOut.write(s[(i + 1) * n : len(s)] + '\n')


fPro = open('/lustre/user/houm/projects/AnnoLnc/RPI_prediction/protein_seq/Homo_sapiens.GRCh38.pep.all.fa')
fFilter = open('/lustre/user/houm/projects/AnnoLnc/RPI_prediction/run4allPros/filtered_all_pro_seq.fa', 'w')

for record in SeqIO.parse(fPro, "fasta"):
	if (not '*' in record.seq) and len(record) > cutoff:
		outputSeq(record, fFilter)

fPro.close()
fFilter.close()