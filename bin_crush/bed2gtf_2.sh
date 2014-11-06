function run_bed2gff {
input_bed=$1
output_gtf=$2
merged_bed=$3
if [ -e $output_gtf ]; then rm -f $output_gtf; fi
if [ ! -e $input_bed ]; then echo "There is no $input_bed"; exit; fi
R --slave <<EOF
options(stringsAsFactors=FALSE)
library("plyr")
bed <- read.table("$input_bed") 

number_exons <- function(x){
#if( length( unique(x[,6] ) ) != 1 ) { 
#	stop( paste( "strand of transcript", unique(x[,4]), " is wrong") ) 
#}
if( unique(x[,6]) %in% c("-", "-1") ) {
	x <- x[order(x[,2], decreasing=TRUE),]
}
if( unique(x[,6]) %in% c("+", "1") ) {
	x <- x[order(x[,2], decreasing=FALSE),]
}
x[,5] <- seq(1, nrow(x) )
x
}
bed <- ddply( bed, .(V4), number_exons )

gtf <- matrix( ".", nrow=nrow(bed) , ncol=9 )
gtf[, 1] <- bed[, 1]
gtf[, 4] <- bed[, 2] 
gtf[, 5] <- bed[, 3]
gtf[, 7] <- bed[, 6]
gtf[, 2] <- rep("bed2gtf", nrow(bed))
gtf[, 3] <- rep("exon", nrow(bed))
gtf[, 6] <- rep(".", nrow(bed))
gtf[, 8] <- rep(".", nrow(bed))

merged_bed <- read.table("$merged_bed", sep="\t")
trans2loci <- vector()
for( i in 1:nrow(merged_bed) ){
	temp.trans <- unlist( strsplit( merged_bed[i, 5], ",") )
	temp.loci <- rep( merged_bed[i, 4], length(temp.trans) )
	names(temp.loci) <- temp.trans
	trans2loci <- c( trans2loci , temp.loci)
}
gtf[, 9] <- paste( 'gene_id "', trans2loci[ bed[, 4] ], '"; transcript_id "', bed[, 4], '"; exon_number"', bed[,5], '";', sep="" )

write.table(gtf, file="$output_gtf", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF
}

run_bed2gff $@
