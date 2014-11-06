#!/usr/bin/env python2.7
import commands

samples = ('SRR070460', 'SRR070448', 'SRR070462')
dSample2Name = {'SRR070448':'FUS', 'SRR070460':'EWSR1', 'SRR070462': 'TAF15'}
fdrs = ('.5', '.1', '.05', '.01', '.001', '.0001')
print 'Number of cross-linking regions\t0.5\t0.1\t0.05\t0.01\t0.001\t0.0001'
for s in samples:
	outLine = [dSample2Name[s]]
	for f in fdrs:
		filePath = s + '_' + f + '/' + s + '_' + f + '_CrossLinkingSites.bed'
		tmp, wcOutput = commands.getstatusoutput('wc -l ' + filePath)
		crosslinkingNum = int(wcOutput.split(' ')[0]) - 1
		# print crosslinkingNum
		outLine.append(str(crosslinkingNum))
	# print outLine
	print '\t'.join(outLine)

