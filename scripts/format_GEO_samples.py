#!/usr/bin/env python2.7

import sys
import re

pBegin = re.compile(r'\d+\. (.*)')

dGSM = {}

for line in sys.stdin.readlines():
	line = line.rstrip('\n')
	
	mBegin = pBegin.match(line)
	# if a new sample begin
	if mBegin:
		# qing kong
		for k in dGSM:
			dGSM[k] = '-'
		dGSM['title'] = mBegin.groups()[0]

	if line.startswith('Source name:'):
		dGSM['source'] = line.partition('\t')[2]
	if line.startswith('Platform:'):
		items = line.split(' ')
		dGSM['platform'] = items[1]
		dGSM['series'] = items[3]
	if line.startswith('Sample'):
		items = line.split()
		dGSM['sample'] = items[2]
		output = '\t'.join([dGSM['sample'], dGSM['series'], dGSM['title'], dGSM['source'], dGSM['platform']])
		print output