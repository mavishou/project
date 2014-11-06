#!/usr/bin/env python2.7

import sys

n = 1
for line in sys.stdin.readlines():	
	line = line.rstrip('\n').split('\t')
	line[3] = line[3] + "_" + str(n)
	print '\t'.join(line)
	n += 1

