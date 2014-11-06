#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
import sys
import modify_gtf as gtf
from collections import defaultdict

ddTrans = defaultdict(int)
curTransId = ''
ok = 1
for line in sys.stdin.readlines():
	li = line.rstrip('\n').split('\t')
	dAnnos = gtf.processAnnotation(li[8])
	lineTransId = dAnnos['transcript_id']
	if lineTransId != curTransId:
		curTransId = lineTransId
		ddTrans[curTransId] += 1

for k, v in ddTrans.items():
	if v > 1:
		sys.stderr.write(k + ' is repeated!\n')
		ok = 0

if ok == 1:
	print 'The file is well sorted by transcripts!'
