#!/usr/bin/python
#from subprocess import call
#call(["Rscript","lib/ZTNB.R","result.filter.rehead.merge","0.05"])
#print "R done"

from lib import getCrosslinking
getCrosslinking.getCrossLinkingMain("2.filter.cluster.bed","2.filter.mutation.bed","2.crosslinking.bed")
print "Done!"

