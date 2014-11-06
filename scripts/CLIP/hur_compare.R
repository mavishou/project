fremoveNo <- function(x){
  return(x[x != '-'])
}

fInterNum <- function(x, y){
  return(length(intersect(x, y)))
}

nm2symbol <- read.table('nm2symbol.txt', sep='\t', stringsAsFactors = F)
tmp <- nm2symbol[, 2]
names(tmp) <- nm2symbol[, 1]
nm2symbol <- tmp
rm(tmp)

sample2ref <- read.table('s2r.txt', sep='\t', stringsAsFactors = F)
sample2ref <- cbind(sample2ref, nm2symbol[sample2ref[, 2]])
sample2ref[, 3] <- as.vector(sample2ref[, 3])
sample2ref[is.na(sample2ref[, 3]), 3] <- '-'

sample2symbol <- sample2ref[, c(1, 3)]
sample2symbol <- unique(sample2symbol)
write.table(sample2symbol, file='HUR_original.txt', sep='\t', quote=F, row.names=F, col.names = F)

aPaper <- sample2symbol[sample2symbol[, 1] == 'HuR_CLIP_A', 2]
bPaper <- sample2symbol[sample2symbol[, 1] == 'HuR_CLIP_B', 2]

aDel <- read.table('GSE28859/SRR189775_deletion_annotation.txt', header = T, stringsAsFactors = F, sep='\t')[, 1]
aDel <- unique(aDel)
aDel  <- fremoveNo(aDel)

aSub <- read.table('GSE28859/SRR189775_substitution_annotation.txt', header = T, stringsAsFactors = F, sep='\t')[, 1]
aSub <- unique(aSub)
aSub <- fremoveNo(aSub)

aIns <- read.table('GSE28859/SRR189775_insertion_annotation.txt', header = T, stringsAsFactors = F, sep='\t')[, 1]
aIns <- unique(aIns)
aIns <- fremoveNo(aIns)


bDel <- read.table('GSE28859/SRR189776_deletion_annotation.txt', header = T, stringsAsFactors = F, sep='\t')[, 1]
bDel <- unique(bDel)
bDel  <- fremoveNo(bDel)

bSub <- read.table('GSE28859/SRR189776_substitution_annotation.txt', header = T, stringsAsFactors = F, sep='\t')[, 1]
bSub <- unique(bSub)
bSub  <- fremoveNo(bSub)

bIns <- read.table('GSE28859/SRR189776_insertion_annotation.txt', header = T, stringsAsFactors = F, sep='\t')[, 1]
bIns <- unique(bIns)
bIns <- fremoveNo(bIns)

aMy <- unique(c(aDel, aIns, aSub))
bMy <- unique(c(bDel, bIns, bSub))


library(VennDiagram)
venn.diagram(list(X=aSub,Y=aIns,Z=aDel),fill=c("red","blue","yellow"),"my_a.tiff")
venn.diagram(list(X=bSub,Y=bIns,Z=bDel),fill=c("red","blue","yellow"),"my_b.tiff")




