# -*- coding: utf-8 -*
#!/usr/bin/python

# 这个脚本用来统计每个位置上碱基的频率，从最后一位开始，向前统计n位
# Author: Hou Mei

import sys

n=20

# result=[{}.fromkeys(['A','C','G','T'],0)] * n
# 初始化result，上面那样初始化只是shadow copy，只产生了一个实例，其他都是指向这个实例的引用，一个变，其他的都变。
# 也可以这样进行初始化：
# result = [fromkeys(['A','C','G','T'],0) for i in range(n)]
# 和下面这种初始化是一样的效果
result=[]
for i in range(n):
	result.append({}.fromkeys(['A','C','G','T'],0))

total=0

for line in sys.stdin.readlines():
	# a sequence or a id?
	if line.startswith('>') == 0:
		# reomove \n
		seq=line.strip()
		seqLen = len(seq)
		total+=1
		i=0
		for p in range(seqLen-n, seqLen):
			result[i][seq[p]]+=1;
			i+=1

print '# Total reads:',total
print 'reverse_position\t#A\t#C\t#G\t#T\t%A\t%C\t%G\t%T'
r=range(0,n)
r.reverse()
for p in r:
	rp=n-p
	dic=result[p]
	num=[dic['A'],dic['C'],dic['G'],dic['T']]
	percent=[round(float(element)/total,2) for element in num]
	output=[rp] + num + percent
	output=[str(el) for el in output]
	print '\t'.join(output)
	
