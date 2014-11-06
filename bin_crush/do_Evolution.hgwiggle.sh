function run_Evolution {
dir=$1
up=$2
down=$3
if [ -d $dir/Evolution.output ]; then rm -fr $dir/Evolution.output ; fi
mkdir $dir/Evolution.output
if [ -e $dir/Evolution.input ]; then rm -fr $dir/Evolution.input ; fi

R --slave <<EOF
options(stringsAsFactors=FALSE)
up <- $up
down <- $down
chr_len <- read.table("blat/hg19.fa.fai", sep="\t", row.names=1)
exons.bed <- read.table("$dir/TransReg.input", sep="\t")
exons.bed[1,3] <- exons.bed[1,3] + 1
strand <- exons.bed[1,6]
chr <- exons.bed[1,1]
if( strand %in% c("+", "1") ){ 
	exons.bed[1,2] <- max( 0, exons.bed[1,2] - up ) 
	exons.bed[1,3] <- min( chr_len[chr, 1], exons.bed[1,3] + down ) } 
if( strand %in% c("-", "0") ) { 
	exons.bed[ 1 ,3] <- min( chr_len[chr, 1] , exons.bed[ 1,3] + up )
	exons.bed[1,2] <- max( 0, exons.bed[1,2] - down ) }
write.table(exons.bed, file="$dir/Evolution.input", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

cd Evolution/phyloP/
../../software_needed/bin/hgWiggle -bedFile=../../$dir/Evolution.input -db=hg19 phyloP46wayAll | grep -P "^[1-9]"  > ../../$dir/Evolution.output/phyloP.All.wig
../../software_needed/bin/hgWiggle -bedFile=../../$dir/Evolution.input -db=hg19 phyloP46wayPlacental | grep -P "^[1-9]" > ../../$dir/Evolution.output/phyloP.Placental.wig
../../software_needed/bin/hgWiggle -bedFile=../../$dir/Evolution.input -db=hg19 phyloP46wayPrimates | grep -P "^[1-9]" > ../../$dir/Evolution.output/phyloP.Primates.wig
cd ../../

R --slave <<EOF
options(stringsAsFactors=FALSE)

region <- read.table("$dir/Evolution.input", sep="\t" )
from <- region[,2]
to <- region[,3]
all.bed <- read.table("$dir/exons.bed", sep="\t" )
colnames(all.bed) <- c("chr", "start", "end", "trans_id", "none", "strand")
IDs <- scan("$dir/SecStructure.output/IDs", what="character")
all.bed <- subset( all.bed, trans_id %in% IDs)
chr <- unique( all.bed[, "chr"] )
strand <- unique( all.bed[, "strand"] )

x <- read.table("blat/hg19.fa.fai", sep="\t" )
chr_length <- x[["V2"]]
names( chr_length ) <- x[["V1"]]

all.scores.list <- list()
score_types <- c("Primates", "Placental", "All")
for( i in  score_types ){
	all.scores.list[[i]] <- read.table( paste("$dir/Evolution.output/phyloP", i, "wig", sep='.'), sep="\t", row.names=1 )
}
DAF <- subset( read.table( bzfile(paste("Evolution/DAF/", chr, ".bz2", sep="") ), sep="\t"), V1 <= to & V1 >= from )
all.scores.list[["DAF"]] <- DAF[,2, drop=FALSE]
rownames( all.scores.list[["DAF"]] ) <- DAF[,1] 
score_types <- c(score_types, "DAF")
write.table(DAF, "$dir/Evolution.output/DAF.wig", sep="\t", row.names=FALSE, col.names=FALSE)

get_mean_score <- function(x, score_type) { mean( all.scores.list[[score_type]][ x, 1 ], na.rm=TRUE ) }

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

        for( i in score_types ) {
                res_score.exon[[t]] <- c( res_score.exon[[t]] , get_mean_score( as.character(temp.exon_pos), i) )
                res_score.up1000[[t]] <- c( res_score.up1000[[t]] , get_mean_score( as.character(temp.up1000), i) )
		nt_Score <- all.scores.list[[ i ]][ as.character(temp.exon_pos), 1 ]
		write( format(nt_Score, digits=3), ncolumns=1, file=paste("$dir/Evolution.output/phyloP", t, i, "nt_score", sep='.') )
        }
}

res <- cbind( t(as.data.frame( res_score.exon ,optional=TRUE ) ), t(as.data.frame( res_score.up1000, optional=TRUE )) )
colnames(res) <- c( paste("exon.", score_types, sep=""), paste("up1000.", score_types, sep="") )
res <- apply( res, c(1,2), function(x){ sprintf("%.2f", x) } )
write.table(res, file="$dir/Evolution.output/Evolution.mean", sep="\t", quote=FALSE, col.names=FALSE)
EOF
}

run_Evolution $@
