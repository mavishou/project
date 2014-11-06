transResult <- read.table('LIN28B.txt', header=T, stringsAsFactors=F, sep = '\t', comment.char = '')

transResult <- transResult[transResult[, 6] >= 2, ]
40742
48733

# 58270
transResult <- transResult[transResult[, 8] / transResult[, 7] >= 3, ]
35818
45989


# sort
aTransResult <- transResult[order(transResult[, 5], decreasing=T), ]
bTransResult <- transResult[order(transResult[, 5], decreasing=T), ]

aGeneIDs <- unique(aTransResult[, 2])
8702
bGeneIDs <- unique(bTransResult[, 2])
10375

> length(intersect(aGeneIDs, bGeneIDs))
[1] 6883

sba <- read.table('sb_LIN28A.txt', stringsAsFactors=F)
sbb <- read.table('sb_LIN28B.txt', stringsAsFactors=F)

ma <- tapply(aTransResult[, 5], as.factor(aTransResult[, 3]), max)
ma <- ma[-1]
mb <- tapply(bTransResult[, 5], as.factor(bTransResult[, 3]), max)
mb <- mb[-1]

sa <- sba[, 4]
names(sa) <- sba[, 1]

sb <- sbb[, 4]
names(sb) <- sbb[, 1]

ia <- intersect(names(sa), names(ma))
cor(sa[ia], ma[ia])

ib <- intersect(names(sb), names(mb))
cor(sb[ib], mb[ib])



