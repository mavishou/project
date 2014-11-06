#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
Author: Mei Hou
This script is used stat the tidy_cluster file. it's the output of tidy_intersect.py
'''

import sys
# import numpy as np
# import pdb
# import cPickle
# import copy
import getopt
from collections import defaultdict

# deal with parameters
bioTypeFile = '/lustre/user/houm/projects/AnnoLnc/bio_type.txt'
tidyFile = ''

opts, args = getopt.getopt(sys.argv[1:], 'b:t:h')
for op, value in opts:
	if op == '-b':
		bioTypeFile = value
	elif op == '-t':
		tidyFile = value
	elif op == '-h':
		print '''Usage: stat_cluster.py -t cluster_tidy_file [-b biotype_file] > output_file

\t-t: the output file tidy_intersect.py
\t-b: the biotype classification file that map a biotype to a protein_coding/pseudogene/lncRNA/sncRNA...'''
		sys.exit()
	else:
		print 'Unknown argument!'
		sys.exit(1)

if tidyFile == '':
	sys.stderr.write('-t must be specified!\n')
	sys.exit(1)

# read in the biotype file
dBioType = {}
bioTypeFile = open(bioTypeFile)
for bt in bioTypeFile.readlines():
	bt = bt.rstrip('\n').split('\t')
	dBioType[bt[0]] = bt[1]


ddFeatureCount = defaultdict(dict)
ddBiotypeCount = defaultdict(dict)
cTrans = {}
cGene = {}
ddChCount = defaultdict(dict)
ddStrandCount = defaultdict(dict)
totalClusters = {}

# read in the cluster tidy file
tidyFile = open(tidyFile)
tidyFile.readline()
for li in tidyFile.readlines():
	li = li.rstrip('\n').split('\t')
	clusterId, ch, strand = li[0:3]
	transId, biotype, feature = li[5:8]
	geneId = li[9]
	
	# check
	if not biotype in dBioType:
		sys.stderr.write(biotype + ' is not in the biotype file!\n')
		sys.exit(1)

	# count total clusters
	totalClusters[clusterId] = 1
	# count chromosome distribution
	ddChCount[ch][clusterId] = 1
	# count strand distribution
	ddStrandCount[strand][clusterId] = 1
	# count transcript num
	cTrans[transId] = 1
	# count gene num
	cGene[geneId] = 1
	# count clusters for every feature
	ddFeatureCount[feature][clusterId] = 1
	# count clusters for every biotype
	ddBiotypeCount[biotype][clusterId] = 1
	

def getCount(mydd):
	'''get count from defaultdict'''
	dOut = {}
	for k, v in mydd.items():
		dOut[k] = len(v)
	return(dOut)

def sortDictAndPrint(myDict, myKey = lambda d: d[1], reverse = True):
	'''input a dict of count, sort it by key, and print it '''
	myDict = sorted(myDict.items(), key = myKey, reverse = reverse)
	for i in myDict:
		print '%s\t%d' % (i[0], i[1])

def dd2Print(mydd, myKey = lambda d: d[1], reverse = True):
	'''conbine getCount() and sortDictAndPrint()'''
	myDict = getCount(mydd)
	sortDictAndPrint(myDict, myKey, reverse)

def toInt(x):
	'''
	input a string. If the string is a putative number, turn it to int, else return itself
	'''
	if x.isdigit():
		return(int(x))
	else:
		return(x)

def printLine():
	print '------------------------------------------------------------------'



cBbioType = getCount(ddBiotypeCount)
cBroadBiotype = defaultdict(int)
for k, v in cBbioType.items():
	cBroadBiotype[dBioType[k]] += v

print '# Clusters that have annotations\t%d' % len(totalClusters)
print '# Transcripts that overlapped with clusters\t%d' % len(cTrans)
print '# Genes that overlapped with clusters\t%d' % len(cGene)
print
printLine()
print
print 'The statistic of clusters for transcript features:'
dd2Print(ddFeatureCount)
print
printLine()
print
print 'The statictic of clusters for biotype groups:'
sortDictAndPrint(cBroadBiotype)
print
printLine()
print
print 'The statictic of clusters for biotypes:'
sortDictAndPrint(cBbioType)
print
printLine()
print
print 'The cluster distribution for both strands: '
dd2Print(ddStrandCount)
print
printLine()
print
print 'The cluster distribution for each chromosomes: '
dd2Print(ddChCount, lambda d: toInt(d[0].split('chr')[1]), False)
print
printLine()



