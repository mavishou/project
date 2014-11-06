#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
Author: Mei Hou
This script is used to count reads for every transcripts, the input is the output of intersectBed result with -a bed -b gtf -wb specified
'''

import sys
import getopt
import bioUtility as bio
from collections import defaultdict

# parameters
collapsed = 0
opts, args = getopt.getopt(sys.argv[1:], 'ch')
for op, value in opts:
	if op == '-c':
		collapsed = 1
	elif op == '-h':
		print 'Usage: countReads4Trans.py [-c] <(standard_in) > output_file\n-c\tThe reads has been collapsed. example: the read name is 1-1000. 1 is its name, 1000 is its duplicated numbers'
		sys.exit(0)

ddCount = defaultdict(int)

if collapsed == 0:
	for line in sys.stdin.readlines():
		li = line.rstrip('\n').split('\t')
		dAnnos = bio.processGTFAnno(li[-1])
		transId = dAnnos['transcript_id']
		ddCount[transId] += 1
else:
	for line in sys.stdin.readlines():
		li = line.rstrip('\n').split('\t')
		readId = li[3]
		dAnnos = bio.processGTFAnno(li[-1])
		readName, count = readId.split('-')
		count = int(count)
		transId = dAnnos['transcript_id']
		ddCount[transId] += count

for k, v in ddCount.items():
	print '%s\t%d' % (k, v)


