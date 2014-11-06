#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
'''
Author: Mei Hou

This script is just for modifying TangX's V4expCaculate.combined.gtf.
It hasn't been tested for other gtf files. 
If you want to use it to deal with other files, make sure the file is sorted by transcript_id

Main Features:  
1. Get introns for every transcripts that have more than 1 exons.
2. Get the 5'UTR and 3'UTR for transcripts that have CDS annotation

Note: 
1. The 9th colunm of the original file is simplified. Only gene_id, transcript_id and exon/intron/CDS_number are retained
2. For some transcripts that have overlaped exons, overlaped exons are conbined when getting introns, but the output of exons are as the original ones.
'''

import sys
import numpy as np
# import pdb
# import cPickle
import copy

class Transcript:
	'''
	geneId
	transID
	exons -> a list of splited gtf exon line lists of transcripts, except for the 9th colunm. Will be sorted by the coordinates
	cdss -> a list of splited gtf CDS line lists of transcripts, except for the 9th colunm
	strand
	generalInfo -> a list of general annotation info of a transcript, contain gtf colunm: 1, 2, 6, 7, 8
	exonCods -> a 2-d numpy array of exon coordinates of a transcripts. sorted
	conbinedExonCods -> a 2-d numpy array of conbined exon coordinates of a transcripts. sorted
	intronCods -> a 2-d numpy array of intron coordinates of a transcripts. sorted
	cdsCods -> a 2-d numpy array of CDS coordinates of a transcripts. sorted
	utr5Cods -> a 2-d numpy array of 5' UTR coordinates of a transcripts. sorted
	utr3Cods -> a 2-d numpy array of 3' UTR coordinates of a transcripts. sorted
	outputList -> a list of output lists of this transcript
	'''
	#pdb.set_trace()
	def __init__(self, geneId, transID, curStrand, exons, cdss):
		self.geneId = geneId
		self.transID = transID
		self.exons = exons
		self.cdss = cdss
		self.strand = curStrand
		self.generalInfo = self.exons[0][0:2] + self.exons[0][5:8]
		self.exons, self.exonCods = self.getCoordiates(self.exons)
		# conbine the overlapped exons
		self.conbinedExonCods = self.getConbineExons()
		self.outputList = []

		self.outputList.extend(self.__getExonOut())
		# self.exons[0].append('test')
		# conbined exon数目大于1才会有intron
		if len(self.conbinedExonCods) > 1:
			self.intronCods = self.getIntronCords()
			self.outputList.extend(self.getIntronOut())

 		# 如果有cds，得到UTR
		if len(self.cdss) > 0:
			self.cdss, self.cdsCods = self.getCoordiates(self.cdss)
			self.utr5Cods, self.utr3Cods = self.getUTR()
			self.outputList.extend(self.__getCDSOut())
			self.outputList.extend(self.getUTROut())
		
		self.getFinalOut()

	def getCoordiates(self, originList):
		'''
		Input: a list of several lines (have been splited to lists) of gtf files. (Only the first 6 colums are requered)
		Return: a 2-d numpy array of sorted coordinates
				the sorted original list
		The sort rule: first by [:, 0] and then by [ :, 1]
		'''
		coordinates = []
		# get the coordinates
		for l in originList:
			coordinates.append([int(s) for s in l[3:5]])
		# turn to np array
		coordinates = np.array(coordinates)
		# get the sort index
		idx = np.lexsort((coordinates[:,1],coordinates[:,0]))
		coordinates = coordinates[idx]
		# orinigal list同时进行排序，否则和坐标不同步了
		# 注意改变这里并不会改变
		# originList = [originList[i] for i in idx]
		originList = copy.deepcopy([originList[i] for i in idx])
		return originList, coordinates

	def conbine2Exons(self, x, y):
		'''
		Input: two overlapped coordinate pairs
		Output: the conbined coordinate pair (in a list)
		Note: the two input pairs must be overlapped.
		'''
		coordinates = np.append(x, y)
		coordinates.sort()
		return [coordinates[0], coordinates[-1]]


	def getConbineExons(self):
		'''
		Main pseudo input: self.exonCods
		Output: self.conbinedExonCods

		'''
		# pdb.set_trace()
		conbinedExonCods = []
		# 用来装那些已经在全局conbine过的，在后面应该把他们去掉
		toBeremoved = []
		for i in range(len(self.exonCods)):
			if not i in toBeremoved:
				toBeremoved.append(i)
				curConbineExon = self.exonCods[i]
				for j in range(i, len(self.exonCods)):
					if not j in toBeremoved:
						if not (curConbineExon[1] < self.exonCods[j][0] or self.exonCods[j][1] < curConbineExon[0]):
							toBeremoved.append(j)
							curConbineExon = self.conbine2Exons(curConbineExon, self.exonCods[j])
				conbinedExonCods.append(curConbineExon)
		conbinedExonCods = np.array(conbinedExonCods)
		return conbinedExonCods

	def getIntronCords(self):
		'''
		Main pseudo input: self.conbinedExonCods
		Output: self.intronCods

		'''		
		intronCods = np.copy(self.conbinedExonCods).reshape(1,-1)[:, 1:-1].reshape(-1,2)
		intronCods[:, 0] = intronCods[:, 0] + 1
		intronCods[:, 1] = intronCods[:, 1] - 1
		return intronCods

	def getUTR(self):
		'''
		Main pseudo input: self.exonCods, self.cdsCods
		Output: self.utr5Cods, self.utr3Cods

		'''	
		utr3 = []
		utr5 = []
		cdscoords = np.copy(self.cdsCods).reshape(1,-1)[0,]
		cdscoords.sort()
		cdsFirst = cdscoords[0]
		cdsLast = cdscoords[-1]
		exoncoords = np.copy(self.exonCods).reshape(1,-1)[0,]
		exoncoords.sort()
		exonFirst = exoncoords[0]
		exonLast = exoncoords[-1]
		if cdsFirst < exonFirst or cdsLast > exonLast:
			sys.stderr.write('The CDS is beyond the exon! ' + self.transID + '\n')
			sys.exit(1)
		# +1和-1不能在检查cds和exon之前，否则，当cds和exon相同时，会报错，而实际上应该是正确的
		cdsFirst -= 1
		cdsLast += 1
		for ec in self.exonCods:
			bf1 = getOverlap(ec, [exonFirst, cdsFirst])
			if bf1 != []:
				utr5.append(bf1)
			bf2 = getOverlap(ec, [cdsLast, exonLast])
			if bf2 != []:
				utr3.append(bf2)
		utr3 = np.array(utr3)
		utr5 = np.array(utr5)
		if self.strand == '-':
			utr5, utr3 = utr3, utr5
		return utr5, utr3


	def getUTROut(self):
		'''
		Main pseudo input: self.utr5Cods, self.utr3Cods, starnd
		Output: the output form of UTR

		'''			
		annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID
		utrOut=[]
		for utr in self.utr5Cods:
			utrOut.append(self.generalInfo[:2] + ["5_UTR"] + [str(s) for s in utr] + self.generalInfo[2:] + [annotation])
		for utr in self.utr3Cods:
			utrOut.append(self.generalInfo[:2] + ["3_UTR"] + [str(s) for s in utr] + self.generalInfo[2:] + [annotation])
		return(utrOut)

	def __getExonOut(self):
		# 如果是负链，exon颠倒
		if self.strand == '-':
			# exonsOut = self.exons[::-1]
			exonsOut = copy.deepcopy(self.exons)[::-1]
		else:
			exonsOut = copy.deepcopy(self.exons)
		for i in range(len(exonsOut)):
			exonNum = i + 1
			annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID + '; exon_number ' + str(exonNum)
			exonsOut[i].append(annotation)
		return(exonsOut)

	def __getCDSOut(self):
		if self.strand == '-':
			cdssOut = copy.deepcopy(self.cdss)[::-1]
		else:
			cdssOut = copy.deepcopy(self.cdss)
		for i in range(len(cdssOut)):
			cdsNum = i + 1
			annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID + '; CDS_number ' + str(cdsNum)
			cdssOut[i].append(annotation)
		return(cdssOut)

	def getIntronOut(self):
		if self.strand == '-':
			myIntronCods = np.copy(self.intronCods)[::-1]
		else:
			myIntronCods = np.copy(self.intronCods)
		intronOut = []
		for i in range(len(myIntronCods)):
			intronNum = i + 1
			annotation = 'gene_id ' + self.geneId + '; transcript_id ' + self.transID + '; intron_number ' + str(intronNum)
			out = self.generalInfo[:2] + ["intron"] + [str(s) for s in myIntronCods[i]] + self.generalInfo[2:] + [annotation]
			intronOut.append(out)
		return(intronOut)
		
	def getFinalOut(self):
		'''
		Main pseudo input: self.outputList
		Ouput: the sorted outputList
		'''
		self.outputList, tmp = self.getCoordiates(self.outputList)
		if self.strand == '-':
			self.outputList = self.outputList[::-1]

def getOverlap(x, y):
	'''
	Input : 2 coordinates pairs
	Return: if they are overlapped, return their overlapped region in a list, otherwise return []
	'''
	if x[1] < y[0] or y[1] < x[0]:
		return []
	else:
		coordinates = np.append(x, y)
		coordinates.sort()
		return coordinates[1:3]

def processAnnotation(annotation):
	'''input: the 9th colume of a gtf file
		output: a dict of annotations 
	'''
	ants = annotation.split(';')
	# 如果以分号结尾，则最后一个元素为空，把空的去掉，并把每个元素开头的空格去掉
	ants = [a.strip() for a in ants if a != '' and a != ' ']
	antsDict = {}
	for a in ants:
		antSplit=a.split(' ')
		##check leng	th
		antsDict[antSplit[0]] = antSplit[1].strip('"')
	return antsDict

def finalOutput(myList):
	'''
		Input: a list output lists
		Action: join lists with '\t' and print them to the standard output
	'''
	for l in myList:
		outputLine = '\t'.join(l)
		print outputLine

####################################################
if __name__ == '__main__':
	curtTransId = ''
	curGeneId = ''
	curStrand = ''
	exons = []
	cdss = []

	# f = open('/lustre/user/houm/projects/AnnoLnc/error.gtf')
	# for line in f.readlines():
	for line in sys.stdin.readlines():
		# pdb.set_trace()
		li = line.rstrip('\n').split('\t')
		# annotation
		# print li
		antsDict = processAnnotation(li[8])
		
		# check whether have gene id and trans id
		if not 'gene_id' in antsDict:
			sys.stderr.write('No gene_id!\n')
			print line
			sys.exit(1)
		if not 'transcript_id' in antsDict:
			sys.stderr.write('No transcript_id!\n')
			print line
			sys.exit(1)

		if antsDict['transcript_id'] != curtTransId and curtTransId !='': ##### 1st
			trans = Transcript(curGeneId, curtTransId, curStrand, exons, cdss)
			finalOutput(trans.outputList)

		if antsDict['transcript_id'] != curtTransId: ##### 2nd
			curtTransId = antsDict['transcript_id']
			curGeneId = antsDict['gene_id']
			curStrand = li[6]
			exons=[]
			cdss=[]
		
		if li[2] == 'exon':
			exons.append(li[:8])
		if li[2] == 'CDS':
			cdss.append(li[:8])

	# the last output
	trans = Transcript(curGeneId, curtTransId, curStrand, exons, cdss)
	finalOutput(trans.outputList)

##test beign
# print str(exons)
# print str(cdss)
# print curStrand
# print curtTransId
# print curGeneId
##test end

##test beign
# w = open('exons.tmp', 'w')
# trans = Transcript(curGeneId, curtTransId, curStrand, exons, cdss)
# # cPickle.dump(trans.exons, w)
# # print len(li)
# # print li[:8]
# # print exons
# print
# print trans.exons
# print
# print trans.exonCods
# print
# print trans.conbinedExonCods
# print
# print trans.intronCods
# print
# print trans.cdss
# print
# print trans.cdsCods
# print
# print trans.utr5Cods
# print
# print trans.utr3Cods
# print
# print trans.outputList
# print
# print trans.transID
# print trans.geneId
# print trans.generalInfo
##test end



