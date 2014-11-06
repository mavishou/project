# =========================function============================
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

# ========================preprocess===========================
# load('exp_36_samples_0715.RData')
# allGenes <- rownames(gfCuff)
# gfCuff <- as.matrix(gfCuff)
# grMy <- as.matrix(grMy)

load('remove_ts_exprs.RData')
genes <- rownames(gfCuff)
# -----------------------filter FPKM------------------------
# co <- 1
# # co <- 1
# coCuff <- apply(gfCuff, 1, function(x) any(x >= co))
# coMy <- apply(grMy, 1, function(x) any(x >= co))
# coAny <- coCuff & coMy
# gfCuff <- gfCuff[coAny, ]
# grMy <- grMy[coAny, ]
# selectGenes <- allGenes[coAny]

# ----------------remove tissu specific genes-----------------
# tsCuff <- fGetTSgenesIndex(gfCuff)
# tsMy <- fGetTSgenesIndex(grMy)
# bothTS <- tsCuff & tsMy
# save(tsCuff, tsMy, file='ts_index.RData')
# gfCuff <- gfCuff[bothTS, ]
# grMy <- grMy[bothTS, ]
# selectGenes <- selectGenes[bothTS]


# ------------pearson correlation--------------
pcCuff <- cor(t(gfCuff))
pcMy <- cor(t(grMy))

# --------------GO preprocess--------------
library(org.Hs.eg.db)
ENSG <- toTable(org.Hs.egENSEMBL)
GO <- toTable(org.Hs.egGO)
GO <- GO[GO[,4] == "BP",]
tmp1Universe <- ENSG[ENSG[, 2] %in% genes, ]
uni4GO <- tmp1Universe[tmp1Universe[, 1] %in% GO[, 1], ]
univ <- unique(uni4GO[, 1])

# =====================GO enrichment================
library(GOstats)
library(stats)

#***
selectGCC <- cor(t(gfCuff))
gcCutOff <- 0.9
#***

genes <- rownames(selectGCC)
selectLnc <- colnames(selectGCC)

ensgList <- entrezList <- lisOfParamObjs <- as.list(1:(2 * ncol(selectGCC)))
ensg2UniqEntrez <- function(ensg, uni4GO){
  unique(uni4GO[uni4GO[, 2] %in% ensg, 1])
}

for(i in 1:ncol(selectGCC)){
  ensgList[[2 * i - 1]] <- genes[selectGCC[, i] >= gcCutOff]
  entrezList[[2 * i - 1]] <- ensg2UniqEntrez(ensgList[[2 * i - 1]], uni4GO)
  lisOfParamObjs[[2 * i - 1]] <- new("GOHyperGParams", geneIds = entrezList[[2 * i - 1]], universeGeneIds = univ, annotation = "org.Hs.eg.db",ontology = "BP")
  ensgList[[2 * i]] <- genes[selectGCC[, i] <= -gcCutOff]
  entrezList[[2 * i]] <- ensg2UniqEntrez(ensgList[[2 * i]], uni4GO)
  lisOfParamObjs[[2 * i]] <- new("GOHyperGParams",geneIds = entrezList[[2 * i]], universeGeneIds = univ, annotation = "org.Hs.eg.db", ontology = "BP")
}

NoEntrezGene <- (rapply(entrezList, length) == 0)
lisOfParamObjs <- lisOfParamObjs[!NoEntrezGene]
library(snowfall)
sfInit(parallel = TRUE, cpus = 30)
hypList <- sfLapply(lisOfParamObjs, hyperGTest)
sfStop()



