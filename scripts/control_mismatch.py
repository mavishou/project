#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
# 用来控制mismatch的数目
# 输入文件是uqniue的reads
import sys
from collections import defaultdict
import getopt
# 从unique_reads脚本调用修正mismatch数目的函数
from unique_reads import fCorrectMismatchNum

# 有一个-m参数，用于指定允许的mismatch数目，默认是0，即不允许mismatch
opts, args = getopt.getopt(sys.argv[1:], 'm:')
allowMismatchNum = 0

for op, value in opts:
	if op == '-m':
		allowMismatchNum = int(value)
	else:
		sys.exit()
# print allowMismatchNum
# print type(allowMismatchNum)

# 统计reads总数
cTotalReads = 0
# 统计每种mismatch数目的defaultdict
ddMismatchStat = defaultdict(int)
# 统计输出了多少reads
cOutputReads = 0

# 如果修正后的mismatch小于等于允许的mismatch，则输出
for line in sys.stdin.readlines():
	cTotalReads += 1
	line = line.rstrip('\n')
	li = line.split('\t')
	correctMismatchNum = fCorrectMismatchNum(li[7])
	# 统计每种mismatch的数目
	ddMismatchStat[correctMismatchNum] += 1
	if correctMismatchNum <= allowMismatchNum:
		print line
		cOutputReads += 1
	

# 统计信息输出到标准错误流
sys.stderr.write('Total reads: ' + str(cTotalReads) + '\n')
sys.stderr.write('Output reads that have equal or less than ' + str(allowMismatchNum) + ' mismatches: ' + str(cOutputReads) + '\n')
sys.stderr.write('Mimatch stat:\n')
for k in sorted(ddMismatchStat):
	sys.stderr.write('# ' + str(k) + ' mismatches: ' + str(ddMismatchStat[k]) + ' reads\n')

