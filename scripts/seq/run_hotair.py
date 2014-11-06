#!/usr/bin/env python2.7
'''This script is used to tidy all protein sequences downloaded from ensembl.
Current filters:
1. sequence no *
2. sequence >50pp'''

from Bio import SeqIO
from math import ceil
import os

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


fPro = open('filtered_all_pro_seq.fa')
fTmp = open('tmp_pro.fa', 'w')
# fResult = open('hotair_all_pro.txt', 'w')
os.system('rm -f hotair_all_pro.txt')

n = 1
cutoff = 1000
for record in SeqIO.parse(fPro, "fasta"):		
		if n > cutoff:
			cmd = 'lncPro HOTAIR_lncRNA.fa tmp_pro.fa >> hotair_all_pro.txt'
			os.system(cmd)
			n = 1
			fTmp.close()
			fTmp = open('tmp_pro.fa', 'w')

		outputSeq(record, fTmp)
		n += 1

fPro.close()
cmd = 'lncPro HOTAIR_lncRNA.fa tmp_pro.fa >> hotair_all_pro.txt'
os.system(cmd)
fTmp.close()

