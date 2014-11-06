# =========================functions------------------------
gcc <- function( x, y ) {
  
  if (!is.vector(x))
    stop("x must be a vector.")
  len = length(x)
  
  if (!is.vector(y) || len != length(y))
    stop("x is a vector; y should be a vector of the same length")
  
  #sort x
  tmp <- sort.int(x, na.last = NA, decreasing = FALSE, index.return = TRUE)
  valuex_x <- tmp$x
  valuey_x <- y[tmp$ix] #y sorted by the rank of x
  
  #sort y
  tmp <- sort.int(y, na.last = NA, decreasing = FALSE, index.return = TRUE)
  valuey_y <- tmp$x
  valuex_y <- x[tmp$ix]  #x sorted by the rank of y
  
  weight <- t(2*seq(1,len,1) - len - 1)
  
  gccxy <- sum(weight*valuex_y)/sum(weight*valuex_x)
  gccyx <- sum(weight*valuey_x)/sum(weight*valuey_y)
  
  #edit by Hou Mei
  #return( data.frame(gccxy, gccyx) )
  gccs <- c(gccxy,gccyx)
  #return the abs max of gccs
  return(max(gccs[abs(gccs)== max(abs(gccs))])) 
  
}

fGetCorGenPer1 <- function(corMatrix1, corMatrix2, cutoff=0.9){
  result <- matrix(0, nrow(corMatrix1), 6)
  rownames(result) <- rownames(corMatrix1)
  colnames(result) <- c('% 1st', '% 2nd', "# Intersect", "# Union", 
                        "# 1st", '# 2nd')
  f1 <- corMatrix1 >= cutoff
  f2 <- corMatrix2 >= cutoff
  f <- f1 + f2
  result[, 3] <- rowSums(f == 2)
  result[, 4] <- rowSums(f >= 1)
  result[, 5] <- rowSums(f1)
  result[, 6] <- rowSums(f2)
  result[, 1] <- result[, 3] / result[, 5]
  result[, 2] <- result[, 3] / result[, 6]
  return(result)
}

fTopOverlap <- function(corMatrix1, corMatrix2, cutoff=100){
  result <- matrix(0, nrow(corMatrix1), 6)
  rownames(result) <- rownames(corMatrix1)
  colnames(result) <- c('% 1st', '% 2nd', "# Intersect", "# Union", 
                        "# 1st", '# 2nd')
  corMatrix1 <- t(apply(corMatrix1, 1, function(x) rank(-x)))
  corMatrix2 <- t(apply(corMatrix2, 1, function(x) rank(-x)))
  f1 <- corMatrix1 <= cutoff
  f2 <- corMatrix2 <= cutoff
  f <- f1 + f2
  result[, 3] <- rowSums(f == 2)
  result[, 4] <- rowSums(f >= 1)
  result[, 5] <- rowSums(f1)
  result[, 6] <- rowSums(f2)
  result[, 1] <- result[, 3] / result[, 5]
  result[, 2] <- result[, 3] / result[, 6]
  return(result) 
}

fGetCorGenPer2 <- function(corMatrix1, corMatrix2, cutoff=-0.9){
  result <- matrix(0, nrow(corMatrix1), 3)
  rownames(result) <- rownames(corMatrix1)
  colnames(result) <- c("# Intersect", "# Union", "Percentage")
  f <- (corMatrix1 <= cutoff) + (corMatrix2 <= cutoff)
  result[, 1] <- rowSums(f == 2)
  result[, 2] <- rowSums(f >= 1)
  result[, 3] <- result[, 1] / result[, 2]
  return(result)
}

fGetDistr <- function(v, intervals=rev(seq(0, 1, 0.05)), type='+'){
  nv <- length(v)
  ni <- length(intervals)
  result <- matrix(0, ni, 1)
  if(type == '+'){
    for(i in 1:ni){
      result[i, 1] <- sum(v >= intervals[i]) / nv
    }
    rownames(result) <- paste('≥', intervals, sep=' ')
  }else{
    for(i in 1:ni){
      result[i, 1] <- sum(v <= intervals[i]) / nv
    }
    rownames(result) <- paste('≤', intervals, sep=' ')
  }
  return(result)
}

getsgene <- function (x, Log = FALSE, Base = 2, AddOne = FALSE, tsThreshold = 0.95, Fraction = TRUE) {
  zeroIndex <- which(apply(x, 1, function(vec) length(which(vec == 0))) == ncol(x))
  if (length(zeroIndex) > 0) {
    x <- x[-zeroIndex, ]
  }
  if (AddOne) {
    x <- x + 1
  }
  if (Log) {
    x <- log(x, Base)
  }
  onets <- function(vec, tsMatrix) {
    tscorematrix <- matrix(0, nrow = dim(tsMatrix)[1], ncol = 2)
    for (i in 1:dim(tsMatrix)[1]) {
      sampleIndex <- which(tsMatrix[i, ] > 0)
      meanvalue <- mean(vec[sampleIndex])
      tscorematrix[i, 1] <- i
      if (meanvalue == 0) {
        tscorematrix[i, 2] <- -100
      }
      else {
        tscorematrix[i, 2] <- 1 - max(vec[-sampleIndex])/meanvalue
      }
    }
    tmax = max(tscorematrix[, 2])
    tmaxidx = tscorematrix[which(tscorematrix[, 2] == tmax), 
                           1][1]
    return(list(tmaxidx = tmaxidx, tmax = tmax))
  }
  tsMatrix <- uniqueTissues(x)
  tt <- apply(x, 1, onets, tsMatrix = tsMatrix)
  tscorematrix <- matrix(0, nrow = dim(x)[1], ncol = 3)
  colnames(tscorematrix) <- c("GeneIndex", "tsmaxscore", "tsmaxidx")
  for (i in 1:dim(x)[1]) {
    tscorematrix[i, 1] <- i
    tscorematrix[i, 2] <- tt[[i]]$tmax
    tscorematrix[i, 3] <- tt[[i]]$tmaxidx
  }
  tscore <- tscorematrix[which(tscorematrix[, 2] >= tsThreshold), 
                         ]
  tsgene <- x[tscore[, 1], ]
  if (Fraction) {
    tsgene <- t(apply(tsgene, 1, function(vec) vec/sum(vec)))
  }
  return(list(tsgene = tsgene, tscore = tscorematrix, uniquets = tsMatrix))
}

uniqueTissues <- function (x) {
  sampleNum <- ncol(x)
  if (is.null(colnames(x))) {
    tsMatrix <- diag(x = 1, nrow = sampleNum, ncol = sampleNum)
  }
  else {
    tsMatrix <- matrix(0, nrow = sampleNum, ncol = sampleNum)
    colnames(tsMatrix) <- colnames(x)
    uniTS <- c("")
    for (i in 1:sampleNum) {
      lastdot <- sapply(gregexpr("\\.", colnames(x)[i]), 
                        tail, 1)
      if (lastdot < 0) {
        curTSName <- colnames(x)[i]
      }
      else {
        if (lastdot > 1) {
          curTSName <- str_sub(colnames(x)[i], 1, lastdot - 
                                 1)
        }
        else {
          curTSName <- colnames(x)[i]
        }
      }
      curTSIndex <- which(uniTS == curTSName)
      if (length(curTSIndex) == 0) {
        if (i == 1) 
          uniTS[1] <- curTSName
        else uniTS <- c(uniTS, curTSName)
        tsMatrix[length(uniTS), i] <- 1
      }
      else {
        tsMatrix[curTSIndex, i] <- 1
      }
    }
    tsMatrix <- tsMatrix[1:length(uniTS), ]
    rownames(tsMatrix) <- uniTS
  }
  tsMatrix
}

fRemoveTSgenes <- function(genes_exprs){
  tsgene <- getsgene(genes_exprs, tsThreshold = 0.95, Fraction = TRUE)
  tsgene_index <- tsgene$tscore[tsgene$tscore[,2]>=0.95,1]
  genes_exprs <- genes_exprs[-tsgene_index, ]
  return(genes_exprs)
}

fGetTSgenesIndex <- function(genes_exprs){
  tsgene <- getsgene(genes_exprs, tsThreshold = 0.95, Fraction = TRUE)
  tsgene_index <- tsgene$tscore[tsgene$tscore[,2]>=0.95,1]
  v <- ! (1:nrow(genes_exprs) %in% tsgene_index)
  return(v)
}

sampleInfo <- read.table('bam/sample_info', stringsAsFactors=F)
samples <- as.vector(t(sampleInfo[, 2:3]))

# ====================cufflinks FPKM=======================
gf4Cuff <- read.table('cuff/cuffnorm2/genes.fpkm_table', header=T, stringsAsFactors=F)
gCuff <- gf4Cuff[, 1]
rownames(gf4Cuff) <- gCuff
gf4Cuff <- gf4Cuff[, -1]

# gf32Cuff <- read.table('cuff/32_samples_cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
# rownames(gf32Cuff) <- gf32Cuff[, 1]
gfCuff <- read.table('cuff/36_classic/genes.fpkm_table', header=T, stringsAsFactors=F)
rownames(gfCuff) <- gfCuff[, 1]
gCuff <- gfCuff[, 1]
# gf32Cuff <- gf32Cuff[, -1]
# gfCuff <- cbind(gf4Cuff, gf32Cuff)
colnames(gfCuff) <- do.call('rbind', strsplit(colnames(gfCuff), '_0'))[, 1]
gfCuff <- gfCuff[, samples]
save(gfCuff, file='gene_fpkm_cuff_classic.RData')

# ===================the count by featurecounts==================
gc32Fc <- read.table('featureCounts/32_sample_gene_reads_count.txt', sep='\t',stringsAsFactors=F, header=T)
rownames(gc32Fc) <- gc32Fc[, 1]
gc32Fc <- gc32Fc[, -(1:6)]
colnames(gc32Fc) <- do.call('rbind', strsplit(colnames(gc32Fc), '[.]'))[, 1]

gc4Fc <- read.table('featureCounts/all_reads_counts.txt', stringsAsFactors = F)
rownames(gc4Fc) <- gc4Fc[, 1]
gc4Fc <- gc4Fc[, -1]
colnames(gc4Fc) <- c('ERR030882', 'ERR030885', 'GSE16256_GSM915328', 'GSE24399_GSM601407')

gcFc <- cbind(gc4Fc, gc32Fc)
gcFc <- gcFc[, samples]
gcFc <- gcFc[gCuff, ]

# total mass
smr <- read.table('featureCounts/32_sample_gene_reads_count.txt.summary', stringsAsFactors=F, header=T)
totalMass32 <- smr[1, -1]
names(totalMass32) <- do.call('rbind', strsplit(names(totalMass32), '[.]'))[, 1]
totalMass32 <- totalMass32[1, ]
totalMass4 <- c(93892894, 101234499, 232204637, 10601172)
names(totalMass4) <- c('ERR030882', 'ERR030885', 'GSE16256_GSM915328', 'GSE24399_GSM601407')
totalMass <- c(totalMass4, totalMass32)
totalMass <- totalMass[samples]

# gene length
gFc <- rownames(gcFc)
geneLength <- read.table('featureCounts/gene_length.txt', stringsAsFactors=F)
geneLength <- geneLength[, -1]
names(geneLength) <- gFc
geneLength <- geneLength[gCuff]

# calculate FPKM
grMy <- gfCuff
grMy[, ] <- 0
for(i in 1:ncol(grMy)){
  grMy[, i] <- gcFc[, i] * 10^9 / (totalMass[i] * geneLength)
}

# save(gcFc, gfCuff, grMy, geneLength, file='exp_36_samples_0714.RData')
save(gcFc, gfCuff, grMy, geneLength, file='exp_36_samples_0715.RData')
save(grMy, file='gene_rpkm_my.RData')

# =================compare correlation=======================
gfCuff <- as.matrix(gfCuff)
grMy <- as.matrix(grMy)

# remove genes that are less than cutoff in any
co <- 0.5
# co <- 1
coCuff <- apply(gfCuff, 1, function(x) any(x >= co))
coMy <- apply(grMy, 1, function(x) any(x >= co))
coAny <- coCuff & coMy
gfCuff <- gfCuff[coAny, ]
grMy <- grMy[coAny, ]
selectGenes <- gCuff[coAny]

# selectGenes <- selectGenes[coAny]

# remove tissu specific genes
tsCuff <- fGetTSgenesIndex(gfCuff)
tsMy <- fGetTSgenesIndex(grMy)
bothTS <- tsCuff & tsMy
save(tsCuff, tsMy, file='ts_index.RData')
gfCuff <- gfCuff[bothTS, ]
grMy <- grMy[bothTS, ]
selectGenes <- selectGenes[bothTS]

save(gfCuff, grMy, file='remove_ts_exprs.RData')

# ----------------compare gene self------------------------
# self pearson correlcation
selfCorrelation <- numeric(length(selectGenes))
names(selfCorrelation) <- selectGenes
for(i in 1:length(selfCorrelation)){
  selfCorrelation[i] <- cor(gfCuff[i, ], grMy[i, ])
}
# selfCorrelation <- selfCorrelation[! is.na(selfCorrelation)]
# write.table(fGetDistr(selfCorrelation), file = 'cor_compare/genes_self_pearson.txt',col.names=F, sep='\t', quote=F)


# remove self pearson less than 0.8
fPearson <- selfCorrelation >= 0.8
gfCuff <- gfCuff[fPearson, ]
grMy <- grMy[fPearson, ]
selectGenes <- rownames(gfCuff)

# self spearman correlation
selfSpearmanCorrelation <- numeric(length(selectGenes))
names(selfSpearmanCorrelation) <- selectGenes
for(i in 1:length(selfSpearmanCorrelation)){
  selfSpearmanCorrelation[i] <- cor(gfCuff[i, ], grMy[i, ], method='spearman')
}
# selfSpearmanCorrelation <- selfSpearmanCorrelation[! is.na(selfSpearmanCorrelation)]
# write.table(fGetDistr(selfSpearmanCorrelation), 
#             file = 'cor_compare/genes_self_spearman.txt',col.names=F, sep='\t', quote=F)

# self gcc cor
selfGcor <- numeric(length(selectGenes))
names(selfGcor) <- selectGenes
for(i in 1:length(selfGcor)){
  selfGcor[i] <- gcc(gfCuff[i, ], grMy[i, ])
}
# write.table(fGetDistr(selfGcor), file = 'cor_compare/genes_self_GCC.txt',col.names=F, sep='\t', quote=F)

selfOut <- cbind(fGetDistr(selfCorrelation), fGetDistr(selfSpearmanCorrelation), fGetDistr(selfGcor))
selfOut <- cbind(rownames(selfOut), selfOut)
colnames(selfOut) <- c('_', 'PCC', "SCC", 'GCC')
write.table(selfOut, file='cor_compare/gene_self_cor_cuff_classic_vs_my.txt', sep='\t', quote=F, row.names=F)

# ----------------calculate correlation-------------------
# pearson correlation
pcCuff <- cor(t(gfCuff))
pcMy <- cor(t(grMy))

# spearman cor
scCuff <- cor(t(gfCuff), method='spearman')
scMy <- cor(t(grMy), method='spearman')

# GCC
gCorCuff <- adjacencymatrix(gfCuff, method="GCC")
gCorMy <- adjacencymatrix(grMy, method="GCC")
diag(gCorCuff) <- 1
diag(gCorMy) <- 1

# ------------------every cor diff compare-------------

# diff of pearson correlation
fDiffCmp <- function(cor1, cor2){
  diff <- abs(cor1 - cor2)
  diff <- diff[upper.tri(diff)]
  result <- fGetDistr(diff, seq(0, 1, 0.05), '-')
  return(result)
}

# diffPc <- abs(pcCuff - pcMy)
# diffPc <- diffPc[upper.tri(diffPc)]
# diffPcDistr <- fGetDistr(diffPc, seq(0, 1, 0.05), '-')
# rm(diffPc)
diffPcDistr <- fDiffCmp(pcMy, pcCuff)
diffScDistr <- fDiffCmp(scMy, scCuff)
diffGcDistr <- fDiffCmp(gCorMy, gCorCuff)

diffDistrOut <- cbind(rownames(diffPcDistr), diffPcDistr, diffScDistr, diffGcDistr)
colnames(diffDistrOut) <- c('-', 'PCC', 'SCC', 'GCC')
write.table(diffDistrOut, file='cor_compare/diff_all_genes.txt', sep='\t', quote=F, row.names=F)

# write.table(fGetDistr(diffPc, seq(0, 1, 0.05), '-'), file = 'diff_all_cors_pearson.txt',col.names=F, sep='\t', quote=F)


summary(diffPc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01890 0.05052 0.08547 0.11260 1.61700 

# spearman
diffSc <- abs(scCuff - scMy)
diffSc <- diffSc[upper.tri(diffSc)]


# ------------compare correlation genes percentage---------------
fOverlapPercentage <- function(cor1, cor2){
  cn <- NULL
#   prefix <- c('0.9', '0.8', '-0.9', '-0.8')
  prefix <- c('0.9', '0.8')
  sufix <- c('# Intersect', '# Union', '%')
  for(i in prefix){
    cn <- cbind(cn, paste(i, sufix, sep=' '))
  } 
  p9 <- fGetCorGenPer1(cor1, cor2)
  p8 <- fGetCorGenPer1(cor1, cor2, 0.8)
#   n9 <- fGetCorGenPer2(cor1, cor2)
#   n8 <- fGetCorGenPer2(cor1, cor2, -0.8)
#   result <- cbind(p9, p8, n9, n8)
  result <- cbind(p9, p8)
  colnames(result) <- cn
  return(result)
}

fOverlapDistr <- function(ol, co=20, type=1){
  if(type == 1){
    tmp <- ol[ol[, 4] >= co, c(1, 2)]
    tmp <- apply(tmp, 1, max)
    distr <- fGetDistr(tmp)
    return(distr)
  }else if(type == 2){
    tmp <- ol[ol[, 4] >= co, c(1, 2)]
    tmp <- c(tmp[, 1], tmp[, 2])
    distr <- fGetDistr(tmp)
    return(distr)
  }else{
    return(fGetDistr(c(ol[, 1], ol[, 2])))
  }
}

fModifyOutput <- function(m, c1 = 'Gene'){
  m <- cbind(rownames(m), m)
  colnames(m)[1] <- c1
  rownames(m) <- NULL
  return(m)
}

# intersect / union percentage, pearson cor >=0.9
# pearson
# pCorPer <- fOverlapPercentage(pcCuff, pcMy)
pCorGenePer9 <- fGetCorGenPer1(pcCuff, pcMy)
# pCorGenePer_9 <- fGetCorGenPer1(pcCuff, pcMy)
# pCorGenePer_8 <- fGetCorGenPer1(pcCuff, pcMy, 0.8)
# pCorGenePer_n9 <- fGetCorGenPer2(pcCuff, pcMy)
# pCorGenePer_n8 <- fGetCorGenPer2(pcCuff, pcMy, -0.8)
# spearman
# sCorPer <- fOverlapPercentage(scCuff, scMy)
# gCorPer <- fOverlapPercentage(gcCuff, gcMy)
sCorGenePer9 <- fGetCorGenPer1(scCuff, scMy)
gCorGenePer9 <- fGetCorGenPer1(gCorCuff, gCorMy)
gCorGenePer95 <- fGetCorGenPer1(gCorCuff, gCorMy, 0.95)

save(pCorGenePer9, sCorGenePer9, gCorGenePer9, gCorGenePer95, file='ol_gene_percent_remove_ts.RData')

fOutputCorGenePer <- function(x, name){
  write.table(fModifyOutput(pCorGenePer9), file=name, sep='\t', quote=F, row.names = F)
}
fOutputCorGenePer(pCorGenePer9, 'cor_compare/ol_compare_pearson_9.txt')
fOutputCorGenePer(sCorGenePer9, 'cor_compare/ol_compare_spearman_9.txt')
fOutputCorGenePer(gCorGenePer9, 'cor_compare/ol_compare_GCC_9.txt')
fOutputCorGenePer(gCorGenePer95, 'cor_compare/ol_compare_GCC_95.txt')

pCorTop100 <- fTopOverlap(pcCuff, pcMy)
pCorTop200 <- fTopOverlap(pcCuff, pcMy, 200)
sCorTop100 <- fTopOverlap(scCuff, scMy)
sCorTop200 <- fTopOverlap(scCuff, scMy, 200)
gCorTop100 <- fTopOverlap(gCorCuff, gCorMy)
gCorTop200 <- fTopOverlap(gCorCuff, gCorMy, 200)
save(pCorTop100, pCorTop200, sCorTop100, sCorTop200, gCorTop100, gCorTop200, file='ol_top_gene_percent_remove_ts')

fOverlapDistr(pCorGenePer9, 50)
fOverlapDistr(pCorGenePer9, type=2)

outOlDisr <- cbind(fOverlapDistr(pCorGenePer9), fOverlapDistr(sCorGenePer9), fOverlapDistr(gCorGenePer9))
outOlDisr <- cbind(rownames(outOlDisr), outOlDisr)
colnames(outOlDisr) <- c('-', 'PCC', 'SCC', 'GCC')
write.table(outOlDisr, file='cor_compare/ol_0.9_compare_1.txt', quote=F, row.names=F, sep='\t')

outOlDisr <- cbind(fOverlapDistr(pCorGenePer9, type=2), fOverlapDistr(sCorGenePer9, type=2), fOverlapDistr(gCorGenePer9, type=2))
outOlDisr <- cbind(rownames(outOlDisr), outOlDisr)
colnames(outOlDisr) <- c('-', 'PCC', 'SCC', 'GCC')
write.table(outOlDisr, file='cor_compare/ol_0.9_compare_2.txt', quote=F, row.names=F, sep='\t')


# fGetDistr(sCorGenePer_9[sCorGenePer_9[, 2] >= 10, 3])
# # GCC
# gCorGenePer_9 <- fGetCorGenPer1(gCorCuff, gCorMy)

# ==============08-11 checking for LiuH=====================
load('exp_36_samples_0715.RData')
liuhGenes <- c('ENSG00000120738', 'ENSG00000134107', 'ENSG00000142408')
liuhGenes %in% rownames(gfCuff)
cor(gfCuff[liuhGenes[1], ], gfCuff[liuhGenes[2], ])
cor(gfCuff[liuhGenes[1], ], gfCuff[liuhGenes[3], ])
cor(grMy[liuhGenes[1], ], grMy[liuhGenes[2], ])
cor(grMy[liuhGenes[1], ], grMy[liuhGenes[3], ])
cor(gfCuff[liuhGenes[1], ], gfCuff[liuhGenes[2], ], method='spearman')
cor(gfCuff[liuhGenes[1], ], gfCuff[liuhGenes[3], ], method='spearman')
  
# get the sample list
sampleInfo <- read.table('sample_info', stringsAsFactors=F, sep='\t')
t <- paste(rep(sampleInfo[, 1], 2), rep(c(1, 2), nrow(sampleInfo)), sep='_')
names(t) <- c(sampleInfo[, 2], sampleInfo[, 3])
exp <- gfCuff[liuhGenes, ]
g <- c('EGR1', 'BHLHE40', 'CACNG8')
rownames(exp) <- g
colnames(exp) <- t[colnames(exp)]
write.table(exp, file='exp.txt', sep='\t', quote=F)
