samples <- scan("clip_samples_GSE.txt", what="character")
starbase <- scan('starbase_GSE.txt', what = 'character')
series <- scan("clip_GEO_series.txt", what = "character")

samples <- unique(samples)
starbase <- unique(starbase)
series <- unique(series)

geo <- unique(c(samples, series))


gsm <- read.table('clip_GEO_samples.txt', sep = "\t", comment.char = '', stringsAsFactors=F)
GSEnotInSamples <- series[! series %in% samples]
add <- matrix("", length(GSEnotInSamples), ncol(gsm))
add[, 2] <- GSEnotInSamples
gsm <- rbind(gsm, add)


infoStarbase <- read.table("starbase_GSE_info.txt", sep = "\t", stringsAsFactors=F)
rbpSB <- tapply(infoStarbase[, 2], as.factor(infoStarbase[, 1]), function(x) {
  y <- unique(x)
  y <- paste(y, collapse='; ')
  return(y)
})
clipTyepSB <- tapply(infoStarbase[, 3], as.factor(infoStarbase[, 1]), function(x) {
  y <- unique(x)
  y <- paste(y, collapse='; ')
  return(y)
})

gsm <- cbind(gsm, '')
gsm <- cbind(gsm, '')
gsm <- cbind(gsm, '')
gsm[, 6] <- rbpSB[gsm[, 2]]
gsm[, 7] <- clipTyepSB[gsm[, 2]]
gsm <- apply(gsm, 2, function(x) {
  x[is.na(x)] <- ''
  return(x)
})
gsm[, 8] <- 0
gsm[gsm[, 6] != '', 8] <- 1

colnames(gsm) <- c('Sample', 'Series', 'Title', 'Source', 'Platform', 'RBP', 'CLIP_type', 'is_in_starBase')
write.table(gsm, file = 'CLIP_GEO_all.txt', sep = '\t', row.names = F, quote = F)

# ---------------human checked------------
checked <- read.table('human_check_clip_GEO.txt', stringsAsFactors=F, header=T, sep = '\t', comment.char = '')
rbps <- checked[, 6]
rbps <- unique(rbps)
rbpSB <- scan('RBPs_starbase.txt', what='character')
rbpOut <- cbind(rbps, as.numeric(rbps %in% rbpSB))
onlyInSB <- rbpSB[! rbpSB %in% rbps]
onlyInSB <- cbind(onlyInSB, 'sb')
rbpOut <- rbind(rbpOut, onlyInSB)


# -----------------
library(limma)
alias2symbolFun <- function(x){
  official <- alias2SymbolTable(x, species = "Hs")
  official[is.na(official)] <- ""
  df <- cbind(x, official)
  df <- as.data.frame(df, stringsAsFactors = F)
  return(df)
}

symbols <- alias2symbolFun(rbpOut[, 1])


# 写出哪些alias map到多个symbol上
alias2MultiSymbolFun <- function(alias, species = "Hs"){
  alias <- as.character(alias)
  species <- "Hs"
  DB <- paste("org", species, "eg", "db", sep = ".")
  ALIAS2EG <- paste("org", species, "egALIAS2EG", sep = ".")
  SYMBOL <- paste("org", species, "egSYMBOL", sep = ".")
  suppressPackageStartupMessages(require(DB, character.only = TRUE))
  isSymbol <- alias %in% Rkeys(get(SYMBOL))
  Symbol <- alias
  Symbol[!isSymbol] <- NA
  OtherAliases <- alias[!isSymbol]
  isAlias <- OtherAliases %in% Rkeys(get(ALIAS2EG))
  #   if (!any(isAlias)) 
  #     return(Symbol)
  OtherAliases <- OtherAliases[isAlias]
  AliasTbl <- toTable(get(ALIAS2EG)[OtherAliases])  
  SymbolTbl <- toTable(get(SYMBOL)[AliasTbl$gene_id])
  # 得到对应到多个gene_id的alias
  dupAlias <- unique(AliasTbl[duplicated(AliasTbl[, 2]), 2])
  t1 <- AliasTbl[AliasTbl[, 2] %in% dupAlias, ]
  t2 <- SymbolTbl[SymbolTbl[, 1] %in% t1[, 1], ]
  m <- merge(t1, t2)
  m <- m[order(m[, 2]), ] # m拿出来挨个检查
  return(m)
}

alias2MultiSymbolFun(rbpOut[, 1])

write.table(symbols, file = 'rbp_official_symbols.txt', sep = '\t', row.names=F, col.names=F, quote=F)
rbpOut <- cbind(rbpOut, symbols[, 2])

write.table(rbpOut, 'rbps.txt', col.names=F, row.names = F, sep = '\t', quote = F)


# ----------GO----------------
rbps <- scan("human_check_RBPs.txt", what = 'character')

GOresult <- read.table('RBPs_GO_terms.txt', sep='\t', header=T, quote='', stringsAsFactors=F)
GOresult <- GOresult[GOresult[, 2] != '', ]
GOresult <- unique(GOresult)
rbps[! rbps %in% GOresult[, 1]]

goTerms <- unique(GOresult[, c(2,3,5)])

genes2GO <- paste(GOresult[, 1], GOresult[, 2], sep = '; ')

out <- matrix(0, nrow(goTerms), length(rbps))
dimnames(out) <- list(goTerms[, 1], rbps)

for(i in 1:length(rbps)){
  tmp <- paste(rbps[i], goTerms[, 1], sep = "; ")
  out[tmp %in% genes2GO, i] <- 1
}
out[out == 0] <- ''

out <- cbind(goTerms[, 2], out)
BPout <- out[goTerms[, 3] == 'biological_process', ]
MFout <- out[goTerms[, 3] == 'molecular_function', ]
CCout <- out[goTerms[, 3] == 'cellular_component', ]

out <- rbind(cbind('BP', BPout), cbind('MF', MFout), cbind('CC', CCout))
out <- cbind(rownames(out), out)
colnames(out)[1:3] <- c("GO accesion", "GO domain", "GO name")

write.table(out, file='RBPs_GO_final.txt', sep = '\t', row.names = F, quote = F)


