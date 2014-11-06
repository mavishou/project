#!/bin_crush/bash
function run_TransReg {
dir=$1 # bed formated
up_dist=$2 # upstream distance from TSS
down_dist=$3 # downstream distance from TES

# checking input and output
if [ -e $dir/TransReg.temp ] ; then rm $dir/TransReg.temp ; fi
if [ -e $dir/TransReg.output ] ; then rm $dir/TransReg.output ; fi
if [ -e $dir/TransReg.input ] ; then rm $dir/TransReg.input ; fi

R --slave <<EOF
options(stringsAsFactors=FALSE)
bed <- read.table("$dir/exons.bed")
TransReg.input <- c(unique(bed[,1]), min(bed[,2]), max(bed[,3]), "query", ".", unique(bed[,6]) )
write.table(t(TransReg.input), file="$dir/TransReg.input", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

bedtools window -sw -l $up_dist -r $down_dist -a $dir/TransReg.input -b TransReg/All_peaks.bed > $dir/TransReg.temp
R --slave <<EOF
options(stringsAsFactors=FALSE)
files <- read.table("TransReg/peak_files.info")
overlaps <- read.table("$dir/TransReg.temp")
gene_region <- read.table("$dir/TransReg.input")
ChIPs <- levels( as.factor( overlaps[, 10] ))

TFBS_table <- matrix( FALSE, nrow=length(ChIPs), ncol=5 )
rownames(TFBS_table) <- ChIPs
colnames(TFBS_table) <- c("up_TSS", "overlap_TSS", "inside_genebody", "overlap_TES", "down_TES")
strand <- gene_region[1,6] 
for( i in 1:nrow(overlaps) ){
	TFBS_start <- overlaps[i, 8]
	TFBS_end <- overlaps[i, 9]
	experiment <- overlaps[i, 10]
	if( strand %in% c("-", "-1") ) { 
		TSS <- gene_region[1, 3] 
		TES <- gene_region[1, 2]
		if ( TSS < TFBS_start ) { TFBS_table[ experiment, "up_TSS"] <- TRUE }
		if ( TES > TFBS_end ) { TFBS_table[ experiment, "down_TES"] <- TRUE }
		if( TES <= TFBS_start & TSS >= TFBS_end ) { TFBS_table[ experiment, "inside_genebody"] <- TRUE }
	}
	else if( strand %in% c("+", "1") ) {
		TSS <- gene_region[1, 2] 
		TES <- gene_region[1, 3]
		if ( TSS > TFBS_end ) { TFBS_table[ experiment, "up_TSS"] <- TRUE }
		if ( TES < TFBS_start ) { TFBS_table[ experiment, "down_TES"] <- TRUE }
		if( TSS <= TFBS_start & TES >= TFBS_end ) { TFBS_table[ experiment, "inside_genebody"] <- TRUE }
	}
	if( TSS >= TFBS_start & TSS <= TFBS_end ) { TFBS_table[ experiment, "overlap_TSS"] <- TRUE }
	if( TES >= TFBS_start & TES <= TFBS_end ) { TFBS_table[ experiment, "overlap_TES"] <- TRUE }
}
output <- cbind( files[ChIPs, ], TFBS_table )
colnames(output) <- c("Cell Line", "Treatment", "Target Protein", "upstream of TSS", "overlap with TSS", "inside genebody", "overlap with TES", "downstream of TES")
write.table( output, file="$dir/TransReg.output", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE )
EOF
sed -i 's/\tFALSE/\t./g' $dir/TransReg.output
}

run_TransReg $@
