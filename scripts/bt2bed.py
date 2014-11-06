#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
Author: Mei Hou
This script is used to change the bowtie format to the bed format
'''

import sys

for line in sys.stdin.readlines():
	li = line.rstrip('\n').split('\t')
	# the begin offset of bt format is 0-based, the same as the bed format
	lOut = li[2:4] + [str(len(li[4]) + int(li[3])), li[0], '.', li[1]]
	print '\t'.join(lOut)

