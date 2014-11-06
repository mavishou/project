#!/usr/bin/env python2.7

from Bio import SeqIO
from collections import defaultdict
import os
from math import ceil

def outputSeq(dFa, id, wOut, n = 70):
	'''a function for output sequence in correct format
	dFa: a dict, key is seq_id, value is sequence(it's str)
	id: the id that want to be output
	wOut: the output file
	'''
	s = dFa[id]
	wOut.write('>' + id + '\n')
	nLine = int(ceil(float(len(s)) / n))
	for i in range(nLine):
		# print s[(i * n) : (i + 1) * n]
		wOut.write(s[(i * n) : (i + 1) * n] + '\n')
	if (i + 1) * n < len(s):
		# print s[(i + 1) * n : len(s)]
		wOut.write(s[(i + 1) * n : len(s)] + '\n')

if __name__ == '__main__':

	fLnc = open('lncRNA_final.fa')
	fPro = open('pro_final.fa')
	fRPI = open('RPI4run.txt')
	# fRPI = open('test.txt')

	dLnc = {}
	dPro = {}

	# read the lncRNA sequence
	for record in SeqIO.parse(fLnc, "fasta"):
		dLnc[record.id] = str(record.seq)

	# read the protein sequence
	for record in SeqIO.parse(fPro, "fasta"):
		dPro[record.id] = str(record.seq)

	dRPI = defaultdict(list)

	# read the RPIs
	for line in fRPI.readlines():
		(lnc, pro) = line.rstrip('\n').split('\t')
		dRPI[lnc].append(pro)

	# prepare and run the lncPro
	for lnc in dRPI:
		# prepare the lncRNA fasta file
		lncFile = lnc + '_lncRNA.fa'
		wLnc = open(lncFile, 'w')
		outputSeq(dLnc, lnc, wLnc)
		wLnc.close()

		# prepare the protein fasta file
		proFile = lnc + '_pro.fa'
		wPro = open(proFile, 'w')
		for pro in dRPI[lnc]:
			outputSeq(dPro, pro, wPro)
		wPro.close()

		# run the lncPro
		resultFile = lnc + '_result.txt'
		cmd = 'lncPro ' + lncFile + ' ' + proFile + ' > ' + resultFile 
		os.system(cmd)



	# after fisnhed lncPro, conbine results
	dResult = defaultdict(defaultdict)
	for lnc in dRPI:
		resultFile = lnc + '_result.txt'
		fResult = open(resultFile)
		for line in fResult.readlines():
			(pro, score) = line.rstrip('\n').split('\t')
			dResult[lnc][pro] = score
		fResult.close()


	# output conbined result
	for lnc in dRPI:
		for pro in dRPI[lnc]:
			out = '\t'.join([lnc, pro, dResult[lnc][pro]])
			print out

	# for k in dRPI:
	# 	print k
	# 	print dRPI[k]
	# print dRPI


	fLnc.close()
	fPro.close()
	fRPI.close()