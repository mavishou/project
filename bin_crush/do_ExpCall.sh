#!/bin_crush/bash

function run_ExpCall {
dir=$1
export PATH=software_needed/cufflinks-2.0.1.Linux_x86_64/:$PATH

if [ -e $dir/ExpCall.all.gtf ]; then rm -f $dir/ExpCall.all.gtf; fi
bash bin_crush/bed2gtf_2.sh $dir/exons.all.bed $dir/ExpCall.all.gtf $dir/merged.bed

if [ -e $dir/ExpCall.output ]; then rm -f $dir/ExpCall.output; fi
if [ -d $dir/ExpCall.temp ]; then rm -rf $dir/ExpCall.temp/; fi

while read i xx
do
	software_needed/cufflinks-2.0.1.Linux_x86_64/cufflinks --no-update-check \
	-o $dir/ExpCall.temp/$i \
	-G $dir/ExpCall.all.gtf \
	ExpCall/bams/$i.bam

done  < ExpCall/cancer/sampleName2ID_131107

while read i len_mean len_sd
do
	software_needed/cufflinks-2.0.1.Linux_x86_64/cufflinks --frag-len-mean $len_mean \
	--frag-len-std-dev $len_sd \
	--no-update-check \
	-o $dir/ExpCall.temp/$i \
	-G $dir/ExpCall.all.gtf \
	ExpCall/bams/$i.bam

done < ExpCall/sample_info

R --slave <<EOF
options(stringsAsFactors=FALSE)
library("plyr")
sample_list <- read.table("ExpCall/sample_list", sep="\t", row.names=1)
raw_FPKM <- list()
loci <- read.table("$dir/merged.bed", sep="\t")[,4]
for( i in c(sample_list[,1], sample_list[,2]) ){
	temp <- read.table( paste("$dir/ExpCall.temp/", i, "/genes.fpkm_tracking", sep=""), row.names=1, sep="\t", header=TRUE )
	raw_FPKM[[i]] <- temp[ loci, "FPKM"]
}
raw_FPKM <- as.data.frame( raw_FPKM )
rownames(raw_FPKM) <- loci
raw_FPKM <- apply( raw_FPKM, c(1,2), as.numeric)
ave_FPKM <- list()
for(i in rownames(sample_list) ){
	ave_FPKM[[i]] <- unlist( alply( raw_FPKM[ , as.vector( as.matrix( sample_list[i, ] ) ) ] , 1, mean ) )
}
ave_FPKM <- as.data.frame(ave_FPKM)
rownames( ave_FPKM ) <- loci

sample_list <- read.table("ExpCall/cancer/sampleName2ID_131107", sep="\t")[,1]
raw_FPKM <- list()
for( i in sample_list ){
	temp <- read.table( paste("$dir/ExpCall.temp/", i, "/genes.fpkm_tracking", sep=""), row.names=1, sep="\t", header=TRUE )
	raw_FPKM[[i]] <- temp[ loci, "FPKM"]
}
raw_FPKM <- as.data.frame( raw_FPKM )
ave_FPKM <- cbind(ave_FPKM, raw_FPKM )

write.table( ave_FPKM, file="$dir/ExpCall.output", sep="\t", quote=FALSE, row.names=TRUE)
EOF

#rm -rf $dir/ExpCall.temp/
}

run_ExpCall $@
