#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
Author: Mei Hou
This script is used to tidy up the result of intersectBed. 
The interSectBed parameters should be seted as: -a xx -b yy -wb [-s], where the -a is a gtf file and -b is a bed6 file
'''

import sys
# import numpy as np
# import pdb
# import cPickle
# import copy
import getopt
# from collections import defaultdict
import modify_gtf as gtf

# deal with parameters
transAnnoFile = '/lustre/user/houm/projects/AnnoLnc/V4_final_transcript_annotation_0401.txt'
clusterAnnoFile = ''
intersectFile=  ''
count = 0
countBiotype = 0
opts, args = getopt.getopt(sys.argv[1:], 't:c:i:h')
for op, value in opts:
	if op == '-t':
		transAnnoFile = value
	elif op == '-c':
		clusterAnnoFile = value
	elif op == '-i':
		intersectFile = value
	elif op == '-h':
		print '''Usage: tidy_intersect.py -g cluster_annotation -i intersect_file [-t transcript_annotaiton] > output_file

\t-c: the output clusters of PARalyzer
\t-i: the output of intersectBed
\t-t: the transcript annotation file abstracted from TangX's V4 gtf file. [default: V4_final_transcript_annotation_0401.txt]'''
		sys.exit()
	else:
		print 'Unknown argument!'
		sys.exit(1)


if clusterAnnoFile == '':
	sys.stderr.write('-c must be specified!\n')
	sys.exit(1)

if intersectFile == '':
	sys.stderr.write('-i must be specified!\n')
	sys.exit(1)

# read in the transcript annotation file
transAnnoFile = open(transAnnoFile)
# remove the first line
transAnnoFile.readline()
dTransAnno = {}
for ta in transAnnoFile.readlines():
	ta = ta.rstrip('\n').split('\t')
	dTransAnno[ta[0]] = ta[1:]

# read in the cluster annotation file
clusterAnnoFile = open(clusterAnnoFile)
# remove the first line
clusterAnnoFile.readline()
dClusterAnno = {}
for ca in clusterAnnoFile.readlines():
	ca = ca.rstrip('\n').split('\t')
	dClusterAnno[ca[4]] = [ca[6]] + ca[9:12]


# tidy the intersect file
print '\t'.join(['ClusterID', 'Chromosome', 'Strand', 'ClusterLength', 'OverlapLength', 'TranscriptID', 'BioType', 'Feature', 'GeneName', 'GeneID', 'ReadCout', 'ConversionLocationCount', 'ConversionEventCount', 'NonConversionEventCount'])
intersectFile = open(intersectFile)
for li in intersectFile.readlines():
	li = li.rstrip('\n').split('\t')
	clusterId = li[12]
	feature = li[2]
	
	# need to count by biotype
	dTransBasicAnno = gtf.processAnnotation(li[8])
	transId = dTransBasicAnno['transcript_id']

	# check
	if not transId in dTransAnno:
		sys.stderr.write(transId + ' is not in the annotation file!\n')
		sys.exit(1)
	# check
	if not clusterId in dClusterAnno:
		sys.stderr.write(clusterId + ' is not in the annotation file!\n')
		sys.exit(1)

	transAnno = dTransAnno[transId]
	biotype = transAnno[1]

	clusterAnno = dClusterAnno[clusterId]
	clusterLen = int(li[11]) - int(li[10])
	intersectLen = int(li[4]) - int(li[3]) + 1
	lOut = [clusterId, li[0], li[6], str(clusterLen), str(intersectLen), transId, biotype, feature, transAnno[2], transAnno[0]] + clusterAnno
	print '\t'.join(lOut)
	




