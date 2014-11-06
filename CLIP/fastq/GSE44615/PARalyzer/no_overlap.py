#!/usr/bin/env python2.7
# -*- coding: utf-8 -*

noOverlapFile = open('/lustre/user/houm/projects/AnnoLnc/CLIP/fastq/GSE44615/PARalyzer/compare_no_overlap_all_groups_LIN28A.txt')
groupFile = open('/lustre/user/houm/projects/AnnoLnc/CLIP/fastq/GSE44615/PARalyzer/LIN28A_groups.txt')

dGroups = {}
for line in noOverlapFile.readlines():
	group = line.split('\t')[3]
	dGroups[group] = 0

print groupFile.readline().rstrip('\n')
for line in groupFile.readlines():
	line = line.rstrip('\n')
	if line.split('\t')[4] in dGroups:
		print line