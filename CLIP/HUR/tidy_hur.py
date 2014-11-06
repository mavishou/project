#!/usr/bin/python

import sys
import re
refseq = []
sample = []
p = re.compile(r'[A-Za-z]')

for line in sys.stdin.readlines():
	line = line.strip()
	# print line
	# line = line.split(' ')
	# print line
	if line.startswith('HuR'):
		tmp = line.split(' ')
		# print tmp
		for t in tmp:
			if t.startswith('HuR'):
				sample.append(t)
	else:
		if p.match(line) and (not line.startswith('Nature')):
		# if p.match(line):
			# print '*************test***************'
			tmp = line.split(' ')
			# print tmp
			for t in tmp:
				if p.match(t):
					refseq.append(t)
# print sample
# print refseq

w1 = open('refseq.txt', 'w')
w2 = open('sample.txt', 'w')

def writeList(l, w):
	for i in range(len(l)):
		w.write(l[i] + '\n')

writeList(refseq, w1)
writeList(sample, w2)
w1.close()
w2.close()