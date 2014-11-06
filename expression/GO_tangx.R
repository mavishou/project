raw_genes_exprs <-read.table("FPKM_table.combine_ES")
gcCutOff <- 0.90
load("/rd/user/houm/v5/grouped_id_list.RData")

# ===============expression filter===============
raw_genes_exprs=as.matrix(raw_genes_exprs)
allGeneType <- c(rep("lncRNA",length(grouped_id_list$lincRNA)),rep("protein_coding",length(grouped_id_list$protein)))
names(allGeneType) <- c(grouped_id_list$lincRNA,grouped_id_list$protein)
f0 <- rownames(raw_genes_exprs) %in% names(allGeneType)
raw_genes_exprs <- raw_genes_exprs[f0,]
gene_type <- allGeneType[rownames(raw_genes_exprs)]
f2=apply(raw_genes_exprs,1,function(x) any(x>0))
genes_exprs <- raw_genes_exprs[f2,]
gene_type <- gene_type[f2]
save(genes_exprs, gene_type, file="130320_exprs_above_0.RData")

# --------------------------
fAbove1 <- apply(genes_exprs, 1, function(x) any(x >= 1))
genes_exprs <- genes_exprs[fAbove1,]
gene_type <- gene_type[fAbove1]

# ===============GCC=====================

library(rsgcc)
gcc_matrix <- adjacencymatrix(genes_exprs, method="GCC")

# ==================GO preprocess======================

library(org.Hs.eg.db)
ENSG <- toTable(org.Hs.egENSEMBL)
GO <- toTable(org.Hs.egGO)
GO <- GO[GO[,4] == "BP",]
tmp1Universe <- ENSG[ENSG[, 2] %in% rownames(genes_exprs), ]
uni4GO <- tmp1Universe[tmp1Universe[, 1] %in% GO[, 1], ]
univ <- unique(uni4GO[, 1])

write.table(univ, file = "univ4GO.txt", row.names = F, col.names = F, quote = F)

# =======================GO enrichment==============
library(GOstats)
library(stats)

fSelect <- gene_type == "lncRNA"
selectGCC <- gcc_matrix[, fSelect]
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
sfInit(parallel = TRUE, cpus = 20)
hypList <- sfLapply(lisOfParamObjs, hyperGTest)
sfStop()

resultList <- lapply(hypList,function(x) summary(x,pvalue=1))
adjpValueList <- lapply(resultList,function(x) p.adjust(x[,2],method="BH"))

enrichCutOff <- 0.01
outPutColNames <- c("lncRNA", "P/N", colnames(resultList[[1]]),"BH")

lncName <- as.character(gl(length(selectLnc), 2, labels = selectLnc))
pn <- rep(c("+","-"),length(selectLnc))
lncName <- lncName[!NoEntrezGene]
pn <- pn[!NoEntrezGene]

#==============output=====================
outPut2 <- NULL
for(i in 1:length(lncName)){
    tmpLen1 <- sum(adjpValueList[[i]] <= 0.01)
	if(tmpLen1 > 0){
		tmp1OutPut <- cbind(lncName[i], pn[i], resultList[[i]][1:tmpLen1, ], adjpValueList[[i]][1:tmpLen1])
		tmp1OutPut <- as.matrix(tmp1OutPut)
		outPut2 <- rbind(outPut2, tmp1OutPut)
	}
}
colnames(outPut2) <- outPutColNames
write.table(outPut2, file = "GCC_.90_merge_ES_all_linc_GO_enrich_result_.01.txt", sep = "\t", row.names = F, quote = F)
save.image("GCC_.90_GO_enrich.RData")
