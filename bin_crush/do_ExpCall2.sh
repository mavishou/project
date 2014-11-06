#!/bin_crush/bash

function run_ExpCall {
dir=$1

if [ ! -e $dir/ExpCall.gtf ] ; then echo "There is no $dir/ExpCall.gtf" ;  exit 1; fi
if [ -e $dir/ExpCall.temp2 ] ; then rm -rf $dir/ExpCall.temp2 ; fi
mkdir $dir/ExpCall.temp2

while read i len_mean len_sd
do
	samtools view -bh ExpCall/bams/$i.bam <(cut -f 1-3 $dir/TransReg.input|sed 's/\t/:/'| sed 's/\t/-/') > $dir/ExpCall.temp2/$i.bam 
	cufflinks --frag-len-mean $len_mean \
	--frag-len-std-dev $len_sd \
	--no-update-check \
	-o $dir/ExpCall.temp2/$i \
	-G $dir/ExpCall.gtf \
	$dir/ExpCall.temp2/$i.bam
	
	lines=`wc -l $dir/ExpCall.temp2/$i/genes.fpkm_tracking | awk '{print $1}'`
	if [ $lines -gt 2 ]; then echo "Gene number greater than 1"; exit 1; fi
done < ExpCall/sample_info


cut -f 1 ExpCall/sample_info | while read i
do
	echo $i 
	samtools view $dir/ExpCall.temp2/$i.bam | cut -f 1 | sort -u | wc -l 
done | paste - - > $dir/ExpCall.temp2/selected_reads_num

R --slave <<EOF
options(stringsAsFactors=FALSE)
sample_list <- read.table("ExpCall/sample_list", sep="\t", row.names=1)
raw_FPKM <- list()
for( i in c(sample_list[,1], sample_list[,2]) ){
	temp <- read.table( paste("$dir/ExpCall.temp2/", i, "/genes.fpkm_tracking", sep=""), sep="\t", header=TRUE )
	raw_FPKM[[i]] <- temp[1, "FPKM"]
}
raw_FPKM <- unlist( raw_FPKM )

all_reads_num <- read.table("ExpCall/mapped_reads", sep="\t", row.names=1 )
selected_reads_num <- read.table("$dir/ExpCall.temp2/selected_reads_num", sep="\t", row.names=1 )
norm_factor <- all_reads_num[ names( raw_FPKM), 1] / selected_reads_num[ names( raw_FPKM), 1]
raw_FPKM <- raw_FPKM / norm_factor

ave_FPKM <- list()
for(i in rownames(sample_list) ){
	ave_FPKM[[i]] <- mean( raw_FPKM[ sample_list[i, ] ] )
}
ave_FPKM <- unlist(ave_FPKM)
write.table(t(ave_FPKM), file="ExpCall.output", sep="\t", quote=FALSE, row.names=FALSE)
EOF

#rm -rf $dir/ExpCall.temp2/
}

run_ExpCall $@
