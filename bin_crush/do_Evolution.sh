#!/bin_crush/bash
function run_Evolution {
dir=$1
if [ -e $dir/Evolution.output ]; then rm -f $dir/Evolution.output ; fi
R --slave <<EOF
options(stringsAsFactors=FALSE)

all.bed <- read.table("$dir/exons.bed", sep="\t" )
colnames(all.bed) <- c("chr", "start", "end", "trans_id", "none", "strand")
chr <- unique( all.bed[, "chr"] )
strand <- unique( all.bed[, "strand"] )

x <- read.table("blat/hg19.fa.fai", sep="\t" )
chr_length <- x[["V2"]]
names( chr_length ) <- x[["V1"]]

score_dirs <- list()
score_dirs[["vertebrate"]] <- "Evolution/phyloP/vertebrate/"
score_dirs[["placentalMammals"]] <- "Evolution/phyloP/placentalMammals/"
score_dirs[["primates"]] <- "Evolution/phyloP/primates/"
score_dirs[["DAF"]] <- "Evolution/DAF/"

all.scores.list <- list()
for( i in names(score_dirs) ) {
	all.scores <- read.table( bzfile(paste(score_dirs[[i]], chr, ".bz2", sep="") ), sep="\t")
	all.scores.list[[i]] <- rep(NA, chr_length[ chr ] )
	all.scores.list[[i]][as.numeric( all.scores[["V1"]] )] <- as.numeric( all.scores[["V2"]] )
}
get_mean_score <- function(x, score_type) { mean( all.scores.list[[score_type]][ as.numeric( x ) ], na.rm=TRUE ) }

res_score.exon <- list()
res_score.up1000 <- list()
for( t in unique( all.bed[, "trans_id"] ) ) {
	sub_bed <- subset( all.bed, trans_id==t )
	temp.exon_pos <- vector()
	for( h in 1:nrow(sub_bed) ){
		temp.exon_pos <- union(temp.exon_pos, seq( sub_bed[h, "start"], sub_bed[h, "end"] ) )
	}

	if( strand %in% c("-", "-1") ) {
		TSS <- max( as.numeric( sub_bed[, "end"] ) )
		temp.up1000 <- seq( TSS, TSS + 1000 )
	}
	if( strand %in% c("+", "1")) {
		TSS <- min( as.numeric( sub_bed[, "start"] ) )
		temp.up1000 <- seq( TSS - 1000, TSS )
	}
	
	for( i in names(score_dirs) ) {
		res_score.exon[[t]] <- c( res_score.exon[[t]] , get_mean_score(temp.exon_pos, i) )
		res_score.up1000[[t]] <- c( res_score.up1000[[t]] , get_mean_score(temp.up1000, i) )
	}
}

res <- cbind( t(as.data.frame( res_score.exon ,optional=TRUE ) ), t(as.data.frame( res_score.up1000, optional=TRUE )) )
colnames(res) <- c( paste("exon.", names(score_dirs), sep=""), paste("up1000.", names(score_dirs), sep="") )
write.table(res, file="$dir/Evolution.output", sep="\t", quote=FALSE)
EOF

}

run_Evolution $@
