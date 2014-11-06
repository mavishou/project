#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
This script is used to get transcript unit info from the tidy file
The input is the tidy_cluster file, the output is the transcript level result
'''

import sys
import numpy as np
# import pdb
# import cPickle
# import copy
import getopt
from collections import defaultdict

# parameters
transCountFile = ''
opts, args = getopt.getopt(sys.argv[1:], 'c:h')
for op, value in opts:
	if op == '-c':
		transCountFile = value
	elif op == '-h':
		print '''

		'''
		sys.exit(0)

if transCountFile == '':
	sys.stderr.write('-c must be specified!\n')
	sys.exit(1)

ddTransCount = defaultdict(int)
transCountFile = open(transCountFile)
for line in transCountFile.readlines():
	li = line.rstrip('\n').split('\t')
	ddTransCount[li[0]] = int(li[1])

def summaryLast():
	global curTransId
	global curBiotype
	global curGeneName
	global curGeneId
	global curCh
	global curStrand 
	global lCurClusterLen
	global lCurOverlaptLen
	global lCurClusterId
	global lCurFeature
	global llCurClusterInfo
	# print llCurClusterInfo
	lCurFeature = np.array(lCurFeature)
	lCurClusterId = np.array(lCurClusterId)
	lCurClusterLen = np.array(lCurClusterLen)
	lCurOverlaptLen = np.array(lCurOverlaptLen)
	llCurClusterInfo = np.array(llCurClusterInfo)
	# print lCurFeatur
	# if exist CDS or UTR, it means that exon is redundant, so remove it
	if 'CDS' in lCurFeature or '3_UTR' in lCurFeature or '5_UTR' in lCurFeature:
		idx = lCurFeature != 'exon'
		lCurClusterId = lCurClusterId[idx]
		llCurClusterInfo = llCurClusterInfo[idx]
		lCurClusterLen = lCurClusterLen[idx]
		lCurOverlaptLen = lCurOverlaptLen[idx]
		lCurFeature = lCurFeature[idx]
	# print idx
	# lCurClusterId = lCurClusterId.tolist()
	# after removed redundant exons, clusters should apear only once in a transcript
	tmp, idx = np.unique(lCurClusterId, return_index = True)
	lCurFeature = tapply(lCurFeature, lCurClusterId, idx)
	lCurOverlaptLen = tapply(lCurOverlaptLen, lCurClusterId, idx)

	lCurClusterId = tmp
	llCurClusterInfo = llCurClusterInfo[idx]
	lCurClusterLen = lCurClusterLen[idx]
	# lCurOverlaptLen = lCurOverlaptLen[idx]

	oClusterId = '; '.join(lCurClusterId)
	oFeature = '; '.join(lCurFeature)
	lCurClusterInfo = [d for d in np.sum(llCurClusterInfo, 0)]
	conversionFrac = str(round(float(lCurClusterInfo[2])/(lCurClusterInfo[3] + lCurClusterInfo[2]), 2))
	oSumClusterInfo = [str(d) for d in lCurClusterInfo]
	tna = np.transpose(llCurClusterInfo)
	oClusterInfo = []
	for a in tna:
		oClusterInfo.append(formatClusterInfo(a))
	# oClusterInfo = np.apply_along_axis(formatClusterInfo, 0, llCurClusterInfo)
	# print oClusterInfo
	# oClusterInfo = np.apply_along_axis(formatClusterInfo, 0, llCurClusterInfo).tolist()
	oClusterLen = '; '.join(lCurClusterLen)
	oOverlapLen = '; '.join(lCurOverlaptLen)
	clusterNum = str(len(lCurClusterId))
	print '\t'.join([curTransId, curGeneId, curGeneName, curCh, curStrand, curBiotype, clusterNum, oClusterId] + 
		oSumClusterInfo + 
		[conversionFrac, str(ddTransCount[curTransId]), oFeature, oClusterLen, oOverlapLen] + 
		oClusterInfo)

def formatClusterInfo(naClusterInfo):
	sl = [str(s) for s in naClusterInfo]
	# print '; '.join(sl)
	return '; '.join(sl)

def tapply(vector, factor, idx):
	'''factor and vector should be np array'''
	result = []
	# tmp, idx = np.unique(factor, return_index = True)
	for i in idx:
		tmp1 = vector[factor == factor[i]]
		result.append(','.join(tmp1))
	return(result)

print '\t'.join(['TranscriptID', 'GeneID', 'GeneName', 'Chromosome', 'Strand', 'BioType', '#Clusters', 'ClusterIDs', 
	'SumReadCout', 'SumConversionLocationCount', 'SumConversionEventCount', 'SumNonConversionEventCount', 'conversionProportion', 'transReadsCount', 'ClusterFeatures', 
	'CLusterLengths', 'OverlapLengths', 'ReadCout', 'ConversionLocationCount', 'ConversionEventCount', 'ConversionEventCount'])

curTransId = ''
curBiotype = ''
curGeneName = ''
curGeneId = ''
curCh = ''
curStrand = ''
lCurClusterLen = []
lCurOverlaptLen = []
lCurClusterId = []
lCurFeature = []
llCurClusterInfo = []

# read file from std in
sys.stdin.readline()
for li in sys.stdin.readlines():
	li = li.rstrip('\n').split('\t')
	clusterId, ch, strand, clusterLen, overlapLen, transId, biotype, feature, geneName, geneId = li[:10]
	lClusterInfo = [int(s) for s in li[10:]]
	# if a cluster is not totally included in a transcript, discard this line
	if clusterLen == overlapLen:
		# summary the last transcript
		if curTransId != '' and transId != curTransId:
			summaryLast()

		# begin a new transcript
		if transId != curTransId:
			curTransId = transId
			curBiotype = biotype
			curGeneName = geneName
			curGeneId = geneId
			curCh = ch
			curStrand = strand
			lCurClusterLen = []
			lCurOverlaptLen = []
			lCurClusterId = []
			lCurFeature = []
			llCurClusterInfo = []

		# every lines should do this
		lCurClusterId.append(clusterId)
		lCurFeature.append(feature)
		llCurClusterInfo.append(lClusterInfo)
		lCurClusterLen.append(clusterLen)
		lCurOverlaptLen.append(overlapLen)

summaryLast()
