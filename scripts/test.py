#!/usr/bin/env python2.7
# -*- coding: utf-8 -*

import sys
from collections import defaultdict
t = 0

ddCount = defaultdict(int)
if t == 1:
	for line in sys.stdin.readlines():
		print line
else:
	for line in sys.stdin.readlines():
		print line + '!'