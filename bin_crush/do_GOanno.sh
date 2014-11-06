function run_GOanno {
dir=$1
if [ -e $dir/GOanno.output.p ] ; then rm -f $dir/GOanno.output.p; fi
if [ -e $dir/GOanno.output.n ] ; then rm -f $dir/GOanno.output.n; fi
R --slave <<EOF
options(stringsAsFactors=FALSE)
library(GOstats)
library(stats)
library(rsgcc)
library("org.Hs.eg.db")
source("bin_crush/gcc4DingY_20121022.R")

gcc_cutoff <- 0.9
adj.p_cutoff <- 0.01

loci <- read.table("$dir/merged.bed", sep="\t")[,4]
FPKM_table <- read.table("GOanno/FPKM_table.univ")
ExpCall.output <- read.table("$dir/ExpCall.output")
GEMatrix <- rbind( ExpCall.output[, colnames(FPKM_table)], FPKM_table )

for( lo in loci){

gcc_res <- getGCCor(lo, as.matrix(GEMatrix) )

uni4GO <- read.table("GOanno/univ4GO.txt") 
univ <- unique(uni4GO[, 1])
ensg2UniqEntrez <- function(ensg, uni4GO){
        unique(uni4GO[uni4GO[, 2] %in% ensg, 1])
}

ensgList <- setdiff( names( gcc_res )[ gcc_res >= gcc_cutoff ], lo )
entrezList <- ensg2UniqEntrez(ensgList, uni4GO)
if( length(entrezList) != 0 ){
	ParamObjs <- new("GOHyperGParams", geneIds = entrezList, universeGeneIds = univ, annotation = "org.Hs.eg.db", ontology = "BP")
	hyper_res <- summary( hyperGTest(ParamObjs), pvalue=1)
	adjpValueList <- p.adjust( hyper_res[,2], method="BH" )
	output.p <- cbind( hyper_res, adjpValueList )[ adjpValueList < adj.p_cutoff, ]
	if( nrow(output.p) > 0 ) {
		colnames( output.p ) <- c( colnames( hyper_res), "BH-adj.Pvalue" )
		write.table( output.p, file=paste("$dir/GOanno.", lo, ".output.p", sep=""), sep="\t", quote=FALSE , row.names=FALSE)
	}
}

ensgList <- setdiff( names( gcc_res )[ gcc_res <= -gcc_cutoff ], lo )
entrezList <- ensg2UniqEntrez(ensgList, uni4GO)
if( length(entrezList) != 0 ){
	ParamObjs <- new("GOHyperGParams", geneIds = entrezList, universeGeneIds = univ, annotation = "org.Hs.eg.db", ontology = "BP")
	hyper_res <- summary( hyperGTest(ParamObjs), pvalue=1)
	adjpValueList <- p.adjust( hyper_res[,2], method="BH" )
	output.n <- cbind( hyper_res, adjpValueList )[ adjpValueList < adj.p_cutoff, ]
	if( nrow(output.n) > 0 ){
		colnames( output.n ) <- c( colnames( hyper_res), "BH-adj.Pvalue" )
		write.table( output.n, file=paste("$dir/GOanno.", lo, ".output.n", sep=""), sep="\t", quote=FALSE , row.names=FALSE)
	}
}
}
EOF
}

run_GOanno $@
