original <- read.table('../original/original_genes.txt', header=T, stringsAsFactors=F)
oGenes <- original[, 1]
oF <- oGenes[original[, 2] > 0]
oE <- oGenes[original[, 3] > 0]
oT <- oGenes[original[, 4] > 0]

setwd('../0.05')
mF <- scan("g_FUS.txt", what="chracters")
mE <- scan("g_EWSR1.txt", what="chracters")
mT <- scan("g_TAF15.txt", what="chracters")
mF <- mF[mF!='-']
mE <- mE[mE!='-']
mT <- mT[mT!='-']

setwd('../original/')
oF <- scan("g_FUS.txt", what="chracters")
oE <- scan("g_EWSR1.txt", what="chracters")
oT <- scan("g_TAF15.txt", what="chracters")
oF <- oF[oF!='-']
oE <- oE[oE!='-']
oT <- oT[oT!='-']

length(intersect(intersect(mF, mE), mT))

length(intersect(oF, mF))
length(intersect(oE, mE))
length(intersect(oT, mT))

library(VennDiagram)
setwd('../0.05')
venn.diagram(list(X=mF,Y=mE,Z=mT),fill=c("red","blue","yellow"),"my.tiff")
setwd('../original/')
venn.diagram(list(X=oF,Y=oE,Z=oT),fill=c("red","blue","yellow"),"original.tiff")

getOverlap <- function(f, e, t){
  all <- unique(c(f, e, t))
  # count table
  ct <- matrix(0, length(all), 3)
  rownames(ct) <- all
  ct[, 1] <- as.numeric(all %in% f)
  ct[, 2] <- as.numeric(all %in% e)
  ct[, 3] <- as.numeric(all %in% t)
  count <- rowSums(ct)
  out <- matrix(0, 3, 2)
  out[1, 1] <- sum(count == 3)
  out[1, 2] <- sum(count > 1)
  out[2, ] <- length(all)
  out[3, ] <- out[1, ] / out[2, ]
  rownames(out) <- c('number', 'total', 'percentage')
  colnames(out) <- c('overlap_in_3', 'overlap_in_at_least_2')
  return(out)
}

get2Overlap <- function(f, e, t){
  all <- unique(c(f, e, t))
  # count table
  ct <- matrix(0, length(all), 3)
  rownames(ct) <- all
  ct[, 1] <- as.numeric(all %in% f)
  ct[, 2] <- as.numeric(all %in% e)
  ct[, 3] <- as.numeric(all %in% t)
  count <- rowSums(ct)
  return(names(count)[count>1])
}

oO <- getOverlap(oF, oE, oT)
oM <- getOverlap(mF, mE, mT)

oOverlap <- get2Overlap(oF, oE, oT)
mOverlap <- get2Overlap(mF, mE, mT)
length(intersect(oOverlap, mOverlap))

rp <- c("SON", "HMGN2", "ERC1", "HNRNPA2B1", "AGPAT6", "NCAPG", "NUDT19", "CACUL1", "ZBTB10", "UCHL5", "TBL1XR1", "DDX17", "RAB21")
all(rp %in% mF)

setwd('../new/')
mF <- scan("g_FUS.txt", what="chracters")
mE <- scan("g_EWSR1.txt", what="chracters")
mT <- scan("g_TAF15.txt", what="chracters")
mF <- mF[mF!='-']
mE <- mE[mE!='-']
mT <- mT[mT!='-']

length(intersect(oF, mF))
length(intersect(oE, mE))
length(intersect(oT, mT))
oM <- getOverlap(mF, mE, mT)
