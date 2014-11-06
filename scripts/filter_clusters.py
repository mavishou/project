#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
This script is used to filter the , the imput is the output of Paralyzer's cluster file
'''
import sys
import getopt

cvfCutoff = 0.25
localtionCutoff = 2
opts, args = getopt.getopt(sys.argv[1:], 'c:l:h')
for op, value in opts:
	if op == '-c':
		cvfCutoff = value
	elif op == '-l':
		localtionCutoff = value
	elif op == '-h':
		print '''Usage: cat INPUT_FILE | filter_clusters.py [-c 0.25 -l 2] > OUTPUT_FILE
\t-c\tThe conversion event proportion cutoff
\t-l\tThe independent location cutoff
		''' 
		sys.exit(0)

firstLine = sys.stdin.readline().rstrip('\n')

totalCluster = 0
leaveCluster = 0

print firstLine
for line in sys.stdin.readlines():
	totalCluster += 1
	line = line.rstrip('\n')
	li = line.split('\t')
	locationCount = int(li[9])
	ConversionCount = float(li[10])
	NonConversionCount = int(li[11])
	conversionFrac = ConversionCount / (ConversionCount + NonConversionCount)
	if locationCount >= localtionCutoff and conversionFrac >= cvfCutoff:
		print line
		leaveCluster += 1


sys.stderr.write('# Total clusters: ' + str(totalCluster) + '\n\n')
sys.stderr.write('Filtered criterion: ' + 'ConversionLocationCount >= ' + str(localtionCutoff) + ', ConversionEventProportion >= ' + str(cvfCutoff) + '\n')
sys.stderr.write('# Left clusters: ' + str(leaveCluster) + '\n')



