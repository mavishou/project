#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
# !!! 输入的bowtie文件一定要先按照名字sort 
# 从bowtie文件中得到unique reads，因为在PAR-CLIP中，T>C和A->G的错配是正常的，要把它们扣除
# 也就是说，一般情况下，strata最高的输出，如果有几个的话就不是unique reads了，这里想换回一下，如果扣掉T>C A>G后，某个reads的strata上来了，那它就是唯一的了

import sys

def fCorrectMismatchNum(mutations):
	'''
	这个函数对mismatch数目进行修正
	输入bowtie文件的最后一列，也就是列出了mutation的那列，扣除正常的mismatch T>C and A>G，返回这个alignment真正mismatch的数目
	如果同时有两种，以次数多的那个为准
	'''
	if mutations == '':
		return 0
	else:
		mutationsLi = mutations.split(',')
		muNum = len(mutationsLi)
		mutationsLi = [m.partition(':')[2] for m in mutationsLi]
		# 从总数里减去正常的mismatch数目
		tToc = mutationsLi.count('T>C')
		aTog = mutationsLi.count('A>G')
		# 如果既有T>C，又有A>G，这是不允许的，减去出现最多的那一个
		muNum = muNum - max(tToc, aTog)
		return muNum

def outputUnique():
	# 前一个read最少的MM数目
	global currentMismathNum
	global uniqueReads
	global currentAlign
	minMismatch = min(currentMismathNum)
	# 如果这个最小的数目是唯一的，那么它是unique read，输出
	countMinMM = currentMismathNum.count(minMismatch)
	if countMinMM == 1:
		idx = currentMismathNum.index(minMismatch)
		print currentAlign[idx]
		uniqueReads += 1

# 如果是main，则是运行此程序，否则是import它，就不需要执行下面的了
if __name__ == '__main__':
	# currentRead是当前正在处理的reads名字
	currentRead = ''
	# currenAlign是当前reads的所有alignment，它是一个list
	currentAlign = []
	# currentMismathNum是当前reads的alignment的mismatch数目
	currentMismathNum = []
	# totalReads是输入的reads总数
	totalReads = 0
	# uniqueReads是unique reads的总数
	uniqueReads = 0

	# 总体思路：从bowtie文件中读取数据，每一行是一个alignment，提取read名字，这时总共有3种情况，可以用currentRead来区别，每种情况需要做的事如下：
	# 1. 第一行 -> currentRead == ''
	# 2. 不是第一行，但是开始了一个新的read: -> currentRead != '' && read_name != currentRead
	# 3. 和上一行read一样 -> read_name == currentRead

	# 需要做的事有如下3点
	# a. 完结前一个reads的各种统计 -> 2
	# b. 开始一个新的read -> 1 2
	# c. 和上一行read一样，read继续 -> 3


	for line in sys.stdin.readlines():
		line = line.rstrip('\n')
		li = line.split('\t')
		# 第7列表示与这行reads相同的其他行的数目，如果它为0，表示肯定是unique reads，直接输出
		if li[6] == 0:
			print line
			uniqueReads += 1
		else:
			if  li[0] != currentRead and currentRead != '': # a
				outputUnique()
			
			if li[0] != currentRead: # b
				totalReads += 1
				currentRead = li[0]
				currentAlign = []
				currentMismathNum = []
				# currentAlign = [line]
				# currentMismathNum = [fCorrectMismatchNum(li[7])]

			currentAlign.append(line)
			currentMismathNum.append(fCorrectMismatchNum(li[7]))

	# 结束之后最后一个还需要输出
	outputUnique()

	# 统计信息输出到标准错误流
	sys.stderr.write('Total reads: ' + str(totalReads) + '\nUnique reads: ' + str(uniqueReads) + '\n')
