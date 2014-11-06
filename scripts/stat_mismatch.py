# -*- coding: utf-8 -*
#!/usr/bin/python
# 用来统计bowtie格式文件mapping，每种错配的比例

import sys
from collections import defaultdict
muTypeDict = defaultdict(int)

for line in sys.stdin.readlines():
	# 前面浙西都是文本处理
	li=line.rstrip('\n').split('\t')
	if li[7] != '':
		mutationsLi = li[7].split(',')
		for m in mutationsLi:
			# print m
			muType = m.partition(':')[2]
			# print muType
			# 这步才是关键
			muTypeDict[muType] += 1
			
# 所有加起来，得到总数，用于算比例
total=sum(muTypeDict.values())
print 'type\tnumber\tpercent'	

# 对key，也就是每种突变类型进行排序
for key in sorted(muTypeDict):
	num = muTypeDict[key]
	# 小心除法中的整数，要转换类型
	outputLi = [key, num, round(float(num)/total, 2)]
	outputLi = [str(e) for e in outputLi]
	print '\t'.join(outputLi)

