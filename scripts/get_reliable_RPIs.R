factor2vector <- function(x){
  for(i in 1:ncol(x)){
    x[, i] <- as.vector(x[, i])
  }
  return(x)
}

library(limma)
diffSymbolFun <- function(x){
  # 输入一串人的gene symbol，输出其被更正或不存在的official name
  official <- alias2SymbolTable(x, species = "Hs")
  official[is.na(official)] <- ""
  ofcIndex <- which(x != official)
  df <- cbind(x, official)
  df <- as.data.frame(df, stringsAsFactors = F)
  df <- df[ofcIndex, ]
  return(df)
}

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

# =======================RAID database====================
library("xlsx")
pr <- read.xlsx2("/Users/Mavis/BaiduYun/AnnoLnc/RAID/RNA-Protein-download.xlsx",1,encoding='UTF-8')

# RAID lncRNA protein interaction
rlnpi <- unique(pr[pr[, 'rna_category'] == 'lncRNA', c(3,7,12,13)])


# ======================guo anyuan==========================
# guo an yuan first
g1 <- read.xlsx2('/Users/Mavis/BaiduYun/AnnoLnc/inter0620-fun.xlsx', 1, encoding='UTF-8')
g1 <- factor2vector(g1)
unique(g1[, 6])
unique(g1[, 'species'])
flag <- (g1[, 6] == "RNA-Protein") & (g1[, 'species'] == 'Homo sapiens') & (g1[, 7] == 'binding')
g1 <- unique(g1[flag, c(4, 5, 2, 8)])

# guo an yuan 2nd -> NPInter
g2 <- read.xlsx2('/Users/Mavis/BaiduYun/AnnoLnc/inter0620-fun.xlsx', 2, encoding='UTF-8')
g2 <- factor2vector(g2)
unique(g2[, 6])
unique(g2[, 'species'])
flag <- (g2[, 6] == "RNA-Protein") & (g2[, 'species'] == 'Homo sapiens') & (g2[, 7] == 'binding')
g2 <- unique(g2[flag, c(4, 5, 2, 8)])

# guo an yuan 3
g3 <- read.xlsx2('/Users/Mavis/BaiduYun/AnnoLnc/inter0620-fun.xlsx', 3, encoding='UTF-8')
g3 <- factor2vector(g3)
unique(g3[, 6])
unique(g3[, 'species'])
flag <- (g3[, 6] == "RNA-Protein") & (g3[, 'species'] == 'Human') & (g3[, 7] == 'binding')
g3 <- unique(g3[flag, c(4, 5, 2, 8)])

# ======================conbine===========================
rlnpi <- cbind(rlnpi, 'RAID')
g1 <- cbind(g1, 'GAY1')
g2 <- cbind(g2, 'GAY2')
g3 <- cbind(g3, 'GAY3')
colnames(rlnpi) <- colnames(g1) <- colnames(g2) <- colnames(g3) <- c('lncRNA', 'Protein', 'Pubmed', 'Description')

all <- rbind(rlnpi, g1, g2, g3)
all <- unique(all)
lnc <- all[, 1]
pro <- all[, 2]

all[, 1] <- toupper(all[, 1])
all[, 2] <- toupper(all[, 2])
rpi <- paste(all[, 1], all[, 2], sep='; ')
pubmed <- tapply(all[, 3], as.factor(rpi), function(x) paste(x, collapse = ' | '))
des <- tapply(all[, 4], as.factor(rpi), function(x) paste(x, collapse = ' ***** '))
src <- tapply(all[, 5], as.factor(rpi), function(x) paste(x, collapse = ' | '))
out <- cbind(pubmed, des, src)
rpi <- do.call('rbind', strsplit(rownames(out), split='; '))
rpi <- as.data.frame(rpi)
out <- cbind(rpi, out)
setwd('../RPI prediction/')
colnames(out) <- c('lncRNA', 'Protein', 'Pubmed', 'Description', 'Source')
write.table(out, file='all4check_rpi.txt', row.names=F, sep='\t', quote=F)

# =========================after check===============================
reliable <- read.table('reliable_rpi_interaction.txt', sep='\t', stringsAsFactors = F)
reliable <- unique(reliable)
lnc <- unique(reliable[, 1])
pro <- unique(reliable[, 2])
write.table(pro, file='pro_final.txt', row.names=F, col.names = F, quote=F)
write.table(reliable, file='RPI_final.txt', row.names=F, col.names = F, quote=F, sep='\t')

# # the disease lncRNA
# dLnc <- scan('/Users/Mavis/BaiduYun/AnnoLnc/RPI prediction/disease_lncRNAs.txt', what='character')
# # lncRNA name in upper case
# dLnc <- toupper(dLnc)

# raid
rLnc <- unique(pr[pr[, 5] == 'lncRNA', c(3, 4)])
rownames(rLnc) <- toupper(rLnc[, 1])
urlLnc <- as.vector(rLnc[, 2])
names(urlLnc) <- toupper(rLnc[, 1])

outLnc <- cbind(lnc, urlLnc[lnc])
write.table(outLnc, file='lncRNAs4check.txt', sep='\t', quote=F, row.names=F, col.names = F)

# protein
alias2MultiSymbolFun(pro)
tmp <- diffSymbolFun(pro)
diff <- tmp[, 2]
names(diff) <- tmp[, 1]
diff[diff == ''] <- 'check'
outPro <- cbind(pro, diff[pro])
outPro[is.na(outPro[, 2]), 2] <- ''
write.table(outPro, file = 'pro4check.txt', sep='\t', quote=F, row.names=F, col.names = F)

# =====================final===============================

