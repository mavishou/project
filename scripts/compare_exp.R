fCompare <- function(x, y, m="pearson"){
  for(i in 2:5){
    print(cor(x[, i], y[, i], method=m))
  }
}

setwd('/lustre/user/houm/projects/AnnoLnc/expression')
# cuff nornal gene count
gcNorm <- read.table('cuff/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
# cuff normal gene fpkm
gfNorm <- read.table('cuff/cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
# gcNorm[, c(2:4)] <- apply(gcNorm[, c(2:4)], 2, as.numeric)
# gfNorm[, c(2:4)] <- apply(gfNorm[, c(2:4)], 2, as.numeric)

# ----------------比较自己的count和fpkm--------------------------
plot(gcNorm[, 2], gfNorm[, 2], pch=20)
for(i in 2:5){
  print(cor(gcNorm[, i], gfNorm[, i], method='spearman'))
}
# 说明同一个组内的count和fpkm其实差不多，也就是说count是对reads长度进行了normalzie的
# 也就是说cuffnorm的作用其实只是在组间进行normalzie
# [1] 0.9772484
# [1] 0.9849432
# [1] 0.9907596
# [1] 0.9897649

# --------------比较normal和mask---------------------------------
# 为了避免组间normalize带来的影响，这里用count来比较
gcMask <- read.table('cuff/mask_rRNA/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
# check whether the order of genes are the same
all(gcMask[, 1] == gcNorm[, 1])
for(i in 2:5){
  print (cor(gcNorm[, i], gcMask[, i]))
}
plot(gcNorm[, 2], gcMask[, 2], pch=20)

# 基本没有区别，说明mask基本上没有影响
# [1] 0.9988739
# [1] 0.9949093
# [1] 0.9999801
# [1] 0.9948179

# Let's check the transcript level
tcNorm <- read.table('cuff/cuffnorm/isoforms.count_table', header=T, stringsAsFactors=F)
tcMask <- read.table('cuff/mask_rRNA/cuffnorm/isoforms.count_table', header=T, stringsAsFactors=F)
# check the order
all(tcNorm[, 1]==tcMask[, 1])
for(i in 2:5){
  print(cor(tcNorm[, i], tcMask[, i]))
}
# transcript level上也是基本上没有区别，以后可以不用mask了
# [1] 0.9983689
# [1] 0.9930249
# [1] 0.9999702
# [1] 0.9917543

# --------------比较normal和correct-------------------------------
### gene count
gcCorrect <- read.table('cuff/correct/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
fCompare(gcNorm, gcCorrect)
# 还是一样的
# [1] 0.9988739
# [1] 0.9949093
# [1] 0.9999801
# [1] 0.9948179

### gene fpkm
gfCorrect <- gcCorrect <- read.table('cuff/correct/cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
fCompare(gfNorm, gfCorrect)
# [1] 0.9982951
# [1] 0.9562654
# [1] 0.9999876
# [1] 0.9528909

### transcript level
tcCorrect <- read.table('cuff/correct/cuffnorm/isoforms.count_table', header=T, stringsAsFactors=F)
fCompare(tcNorm, tcCorrect)
# [1] 0.9983689
# [1] 0.9930249
# [1] 0.9999702
# [1] 0.9917543

# --------------比较normal和no_length-------------------------------
### gene count
gcNl <- read.table('cuff/no_length/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
fCompare(gcNorm, gcNl)
# [1] 0.9988739
# [1] 0.9949093
# [1] 0.9999801
# [1] 0.9948179
plot(gcNorm[, 2], gcNl[, 2], pch=20)

### gene fpkm
gfNl <- read.table('cuff/no_length/cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
# check the order
all(gfNl[, 1] == gfNorm[, 1])
fCompare(gfNorm, gfNl)
# [1] 0.1080504
# [1] 0.4147376
# [1] 0.1051083
# [1] 0.6595122
fCompare(gfNorm, gfNl, 'spearman')
# [1] 0.9767756
# [1] 0.9843196
# [1] 0.9904253
# [1] 0.9894595

plot(gfNorm[, 2], gfNl[, 2], pch=20)

# ---------------比较norm和HTSeq-----------------------------
### gene count
gcHtseq <- read.table('HTSeq/gene.txt', sep='\t', stringsAsFactors=F)
# check ID
gCuff <- gcNorm[, 1]
gHtseq <- gcHtseq[, 1]
all(gCuff %in% gHtseq)
# setdiff(gHtseq, gCuff)
gcHtseq <- gcHtseq[gcHtseq[, 1] %in% gCuff, ]
# check order
all(gcHtseq[, 1] == gcNorm[, 1])

fCompare(gcNorm, gcHtseq)
# [1] 0.6106263
# [1] 0.3612847
# [1] 0.3582761
# [1] 0.9885178
fCompare(gcNorm, gcHtseq, 'spearman')
# [1] 0.9195368
# [1] 0.9267516
# [1] 0.9134145
# [1] 0.9461353

plot(gcNorm[, 2], gcHtseq[, 2], pch=20, xlim=c(0, 8e5), ylim=c(0, 8e5), 
     xlab="Cuff", ylab="HTSeq", main="ERR030882 genes")
text(6e5, 6e5, 'Pearson cor: 0.61', cex=1.2)
text(6e5, 7e5, 'Spearman cor: 0.92', cex=1.2)

# 找出在cuff中为0，却在HTSeq中reads数目很多的
abNormal <- gCuff[gcNorm[, 2] == 0 & gcHtseq[, 2] > 100000]
rownames(gcNorm) <- rownames(gcNorm2) <- rownames(gcHtseq) <- gCuff
gcNorm[abNormal, ]

plot(gcNorm[, 5], gcHtseq[, 5], pch=20, 
     xlim=c(0, 4e5), ylim=c(0, 4e5), 
     xlab="Cuff", ylab="HTSeq", main="GSE24399_GSM601407 genes")
text(2e5, 2e5, 'Pearson cor: 0.99', cex=1.2)
text(2e5, 3e5, 'Spearman cor: 0.95', cex=1.2)

### transcript count
tcHtseq <- read.table('HTSeq/transcript.txt', sep='\t', stringsAsFactors=F)
tcHtseq <- tcHtseq[tcHtseq[, 1] %in% tcNorm[, 1], ]
# check ordr
all(tcNorm[, 1] == tcHtseq[, 1])

fCompare(tcNorm, tcHtseq)
# [1] 0.1980898
# [1] 0.1412311
# [1] 0.09203862
# [1] 0.6937974

fCompare(tcNorm, tcHtseq, 'spearman')
# [1] 0.5125536
# [1] 0.5286943
# [1] 0.5552942
# [1] 0.4875728

plot(tcNorm[, 2], tcHtseq[, 2], pch=20, 
     xlim=c(0, 8e5), ylim=c(0, 8e5), 
     xlab="Cuff", ylab="HTSeq", main="ERR030882 transcripts")
text(4e5, 4e5, 'Pearson cor: 0.20', cex=1.2)
text(4e5, 5e5, 'Spearman cor: 0.51', cex=1.2)

# ---------------------------Norm2--------------------------------
gcNorm2 <- read.table('cuff/cuffnorm2/genes.count_table', header=T, stringsAsFactors=F)
gfNorm2 <- read.table('cuff/cuffnorm2/genes.fpkm_table', header=T, stringsAsFactors=F)

# -----------------------featureCount---------------------------
gcFc <- read.table('featureCounts/all_gene_counts.txt', sep='\t', stringsAsFactors=F)
# genes of featurecount
gFc <- gcFc[, 1]
all(gCuff %in% gFc)
rownames(gcFc) <- gFc
gcFc <- gcFc[, ]
gcFc <- gcFc[gCuff, ]
# check order
all(gcFc[, 1] == gcNorm[, 1])

# compare featureCount with HTSeq
fCompare(gcHtseq, gcFc)
> fCompare(gcHtseq, gcFc)
# [1] 0.9930448
# [1] 0.9989145
# [1] 0.996702
# [1] 1

# calculate FPKM using reads count
# total mass
totalMass <- c(4.92436e+07, 5.08075e+07, 1.15251e+08, 1.20944e+07)
# gene length
geneLength <- read.table('featureCounts/gene_length.txt', stringsAsFactors=F)
geneLength <- geneLength[, -1]
names(geneLength) <- gFc
geneLength <- geneLength[gCuff]

gfMy <- gfNorm
gfMy[, 2:5] <- 0
for(i in 2:5){
  gfMy[, i] <- gcFc[, i] * 10^9 / (totalMass[i-1] * geneLength)
}

fCompare(gfMy, gfNorm2)
# [1] 0.02718042
# [1] 0.1237496
# [1] 0.01356684
# [1] 0.7207268

fCompare(gfMy, gfNorm2, 'spearman')
# [1] 0.8766704
# [1] 0.8983167
# [1] 0.8718634
# [1] 0.9215698

fCompare(gfMy, gfNorm)
# [1] 0.02718064
# [1] 0.1237496
# [1] 0.01356685
# [1] 0.7207269

fCompare(gfMy, gfNorm, 'spearman')
# [1] 0.8766291
# [1] 0.8982885
# [1] 0.8718319
# [1] 0.9215441

plot(gfNorm2[, 2], gfMy[, 2], pch=20)
plot(gfNorm2[, 2], gfMy[, 2], pch=20, xlim=c(0, 1000), ylim=c(0, 1000))

# ------------compare fc and HTSeq and Norm2 ----------------------
samples <- c('ERR030882', 'ERR030885', 'GSE16256_GSM915328', 'GSE24399_GSM601407')
fCompare(gcFc, gcHtseq)
# [1] 0.9930448
# [1] 0.9989145
# [1] 0.996702
# [1] 1
png(file="1.png", height = 500, width = 500, pointsize = 15)
par(mfrow=c(2,2))
for(i in 2:5){
  plot(gcFc[, i], gcHtseq[, i], pch=20, xlab='FeatureCounts', ylab='HTSeq', 
       main=paste(samples[i-1], 'fragments count', sep=' '))
}

for(i in 2:5){
  plot(gcNorm2[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm (not normalized for library size)', 
       main=paste(samples[i-1], 'fragments count', sep=' '))
}

i=2
plot(gcNorm2[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm (not normalized for library size)', 
     xlim=c(0, 6e5), ylim=c(0, 6e5), main=paste(samples[i-1], 'fragments count', sep=' '))
i=3
plot(gcNorm2[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm (not normalized for library size)', 
     xlim=c(0, 6e5), ylim=c(0, 6e5), main=paste(samples[i-1], 'fragments count', sep=' '))
i=4
plot(gcNorm2[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm (not normalized for library size)', 
     xlim=c(0, 150000), ylim=c(0, 150000), main=paste(samples[i-1], 'fragments count', sep=' '))
i=5
plot(gcNorm2[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm (not normalized for library size)', 
     xlim=c(0, 1e5), ylim=c(0, 1e5), main=paste(samples[i-1], 'fragments count', sep=' '))

for(i in 2:5){
  plot(gcNorm[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm', 
       main=paste(samples[i-1], 'fragments count', sep=' '))
}

i=2
plot(gcNorm[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm', 
     xlim=c(0, 6e5), ylim=c(0, 6e5), main=paste(samples[i-1], 'fragments count', sep=' '))
i=3
plot(gcNorm[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm', 
     xlim=c(0, 6e5), ylim=c(0, 6e5), main=paste(samples[i-1], 'fragments count', sep=' '))
i=4
plot(gcNorm[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm', 
     xlim=c(0, 150000), ylim=c(0, 150000), main=paste(samples[i-1], 'fragments count', sep=' '))
i=5
plot(gcNorm[, i], gcFc[, i], pch=20, ylab='FeatureCounts', xlab='CuffNorm', 
     xlim=c(0, 1e5), ylim=c(0, 1e5), main=paste(samples[i-1], 'fragments count', sep=' '))

# ------------------compare FPKM and counts---------------------
i=2
plot(gcNorm2[, i], gfNorm2[, i], pch=20, ylim=c(0, 1000), xlim=c(0, 1e5), 
     xlab='Fragments count', ylab='FPKM', main=paste(samples[i-1], 'Cuffnorm', sep=' '))
i=3
plot(gcNorm2[, i], gfNorm2[, i], pch=20, xlim = c(0, 1e5), ylim = c(0, 5000),
     xlab='Fragments count', ylab='FPKM', main=paste(samples[i-1], 'Cuffnorm', sep=' '))

i=4
plot(gcNorm2[, i], gfNorm2[, i], pch=20, ylim = c(0, 1e4), xlim=c(0, 2e5),
     xlab='Fragments count', ylab='FPKM', main=paste(samples[i-1], 'Cuffnorm', sep=' '))

i=5
plot(gcNorm2[, i], gfNorm2[, i], pch=20, 
     xlab='Fragments count', ylab='FPKM', main=paste(samples[i-1], 'Cuffnorm', sep=' '))

# -------------compare my FPKM and cuff FPKM--------------

i=2
plot(gfNorm2[, i], gfMy[, i], pch=20, xlab='CuffNorm (not normalized for library size)', ylab='HM', 
     xlim=c(0, 1500), ylim=c(0, 1500),
     main=paste(samples[i-1], 'FPKM', sep=' '))

i=3
plot(gfNorm2[, i], gfMy[, i], pch=20, xlab='CuffNorm (not normalized for library size)', ylab='HM', 
     xlim=c(0, 2500), ylim=c(0, 2500),
     main=paste(samples[i-1], 'FPKM', sep=' '))

i=4
plot(gfNorm2[, i], gfMy[, i], pch=20, xlab='CuffNorm (not normalized for library size)', ylab='HM', 
     xlim=c(0, 2500), ylim=c(0, 2500),
     main=paste(samples[i-1], 'FPKM', sep=' '))

i=5
plot(gfNorm2[, i], gfMy[, i], pch=20, xlab='CuffNorm (not normalized for library size)', ylab='HM', 
     xlim=c(0, 3000), ylim=c(0, 3000),
     main=paste(samples[i-1], 'FPKM', sep=' '))


for(i in 2:5){
  plot(gfNorm2[, i], gfMy[, i], pch=20, xlab='CuffNorm (not normalized for library size)', ylab='HM', 
       xlim=c(0, 500), ylim=c(0, 500),
       main=paste(samples[i-1], 'FPKM', sep=' '))
}


# -----------------------cuffdiff output-----------------------------
gcCuffRaw <- read.table('norm3_raw_counts.txt', stringsAsFactors = F)
all(gcCuffRaw[, 1] == gcNorm2[, 1])
fCompare(gcFc, gcCuffRaw)
# [1] 0.6002613
# [1] 0.3629053
# [1] 0.3543792
# [1] 0.9885178

fCompare(gcFc, gcCuffRaw, 'spearman')
# [1] 0.9252406
# [1] 0.9320169
# [1] 0.9172879
# [1] 0.9461658


par(mfrow=c(2,2))
for(i in 2:5){
  plot(gcFc[, i], gcCuffRaw[, i], pch=20, xlab='FeatureCounts', ylab='Cufflinks raw', 
       main=paste(samples[i-1], 'fragments count', sep=' '))
}
i=2
plot(gcFc[, i], gcCuffRaw[, i], pch=20, xlab='FeatureCounts', ylab='Cufflinks raw', 
     main=paste(samples[i-1], 'fragments count', sep=' '),
     xlim=c(0, 6e5), ylim=c(0, 6e5))

i=3
plot(gcFc[, i], gcCuffRaw[, i], pch=20, xlab='FeatureCounts', ylab='Cufflinks raw', 
     main=paste(samples[i-1], 'fragments count', sep=' '),
     xlim=c(0, 6e5), ylim=c(0, 6e5))

i=4
plot(gcFc[, i], gcCuffRaw[, i], pch=20, xlab='FeatureCounts', ylab='Cufflinks raw', 
     main=paste(samples[i-1], 'fragments count', sep=' '),
     xlim=c(0, 6e5), ylim=c(0, 6e5))

i=5
plot(gcFc[, i], gcCuffRaw[, i], pch=20, xlab='FeatureCounts', ylab='Cufflinks raw', 
     main=paste(samples[i-1], 'fragments count', sep=' '), 
     xlim=c(0, 8e4), ylim=c(0, 8e4))


# -----------not correct for effective length-------------------
gcEffective <- read.table('genes.count_table', header = T, stringsAsFactors = F)
gfEffective <- read.table('genes.fpkm_table', header = T, stringsAsFactors = F)

par(mfrow=c(2,2))
for(i in 2:5){
  plot(gfEffective[, i], gfMy[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='HM',
       main=paste(samples[i-1], 'FPKM', sep=' '))
}

i=2
plot(gfEffective[, i], gfMy[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='HM', 
     xlim=c(0, 1500), ylim=c(0, 1500),
     main=paste(samples[i-1], 'FPKM', sep=' '))

i=3
plot(gfEffective[, i], gfMy[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='HM', 
     xlim=c(0, 2500), ylim=c(0, 2500),
     main=paste(samples[i-1], 'FPKM', sep=' '))

i=4
plot(gfEffective[, i], gfMy[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='HM', 
     xlim=c(0, 2500), ylim=c(0, 2500),
     main=paste(samples[i-1], 'FPKM', sep=' '))

i=5
plot(gfEffective[, i], gfMy[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='HM', 
     xlim=c(0, 3000), ylim=c(0, 3000),
     main=paste(samples[i-1], 'FPKM', sep=' '))

for(i in 2:5){
  plot(gfEffective[, i], gfMy[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='HM', 
       xlim=c(0, 500), ylim=c(0, 500),
       main=paste(samples[i-1], 'FPKM', sep=' '))
}


# ------------edgeR--------------------------
dge <- DGEList(counts=gcFc[, 2:5], lib.size=totalMass)
# dge$genes$length <- geneLength
gfEdgeR <- rpkm(dge, gene.length=geneLength)



# ---------------RPKM-----------------------
gcrFc <- read.table('all_reads_counts.txt', stringsAsFactors = F)

rownames(gcrFc) <- gFc
gcrFc <- gcrFc[gCuff, ]
# check order
all(gcrFc[, 1] == gcNorm[, 1])

fCompare(gcrFc, gcFc)
# [1] 0.9980141
# [1] 0.9992908
# [1] 0.9977601
# [1] 1

# calculate FPKM using reads count
# total mass
totalMass2 <- c(134740145, 143922366, 363856013, 17407482)
totalMass3 <- c(93892894, 101234499, 232204637, 10601172)

grMy <- gfNorm
grMy[, 2:5] <- 0
for(i in 2:5){
  grMy[, i] <- gcrFc[, i] * 10^9 / (totalMass3[i-1] * geneLength)
}

dge <- DGEList(counts=gcrFc[, 2:5], lib.size=totalMass3)
grEdgeR <- rpkm(dge, gene.length=geneLength)


# ---------------effective length and no effective length-----------------
for(i in 2:5){
  plot(gfEffective[, i], gfNorm2[, i], pch=20, xlab='Cuffquant (not correlated for effective length)', ylab='Cuffquant', 
       xlim=c(0, 3000), ylim=c(0, 3000),
       main=paste(samples[i-1], 'FPKM', sep=' '))
}




