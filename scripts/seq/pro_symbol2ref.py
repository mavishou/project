#!/usr/bin/env python2.7

from collections import defaultdict

fPro = open('pro_final.txt')
fBiomart = open('biomart_pro_symbol2ref.txt')

lPro = []

for line in fPro.readlines():
	line = line.rstrip('\n')
	lPro.append(line)

dSymbol2ref = defaultdict(str)

for line in fBiomart.readlines():
	line = line.rstrip('\n').split('\t')
	if line[1] != '':
		if not line[0] in dSymbol2ref:
			dSymbol2ref[line[0]] = line[1]

for pro in lPro:
	out = '\t'.join([pro, dSymbol2ref[pro]])
	print out