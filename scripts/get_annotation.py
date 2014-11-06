#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
Author: Mei Hou
Get the transcript annotation from TangX's V4expCaculate.combined.gtf
Include: 'transcript_id', 'gene_id', 'gene_biotype', 'transcript_type', 'source', 'gene_name', 'transcript_name'
Note: if a annotation type is not existed for a transcript, it is '-'
'''

import sys
import modify_gtf as gtf

curTransId = ''
tAnnoType = ('transcript_id', 'gene_id', 'gene_biotype', 'transcript_type', 'source', 'gene_name', 'transcript_name')
print '\t'.join(tAnnoType)
for line in sys.stdin.readlines():
	li = line.rstrip('\n').split('\t')
	# 只看exon的，以防其他UTR这些注释不全
	if li[2] == 'exon':
		dAnnos = gtf.processAnnotation(li[8])
		dAnnos['source'] = li[1]
		if curTransId != dAnnos['transcript_id']:
			curTransId = dAnnos['transcript_id']
			lOut = []
			for at in tAnnoType:				
				lOut.append(dAnnos.setdefault(at, '-'))
			print '\t'.join(lOut)



