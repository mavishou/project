function run_disease {
dir=$1
up_dist=$2 # upstream distance from TSS
down_dist=$3 # downstream distance from TES

if [ -e $dir/Disease.input ]; then rm -f $dir/Disease.input ; fi
if [ -e $dir/Disease.temp1 ] ; then rm -f $dir/Disease.temp1 ; fi
if [ -e $dir/Disease.temp2 ] ; then rm -f $dir/Disease.temp2 ; fi
if [ -e $dir/Disease.output ] ; then rm -f $dir/Disease.output ; fi
R --slave <<EOF
options(stringsAsFactors=FALSE)
bed <- read.table("$dir/exons.bed")
strand <- unique( bed[,6] )
if( strand %in% c("+", "1")) {
	Disease.input <- c(unique(bed[,1]), min(bed[,2]) - $up_dist, max(bed[,3]) + $down_dist, "query", ".", unique(bed[,6]) )
}
if( strand %in% c("-", "-1")) {
	Disease.input <- c(unique(bed[,1]), min(bed[,2]) - $down_dist, max(bed[,3]) + $up_dist, "query", ".", unique(bed[,6]) )
}
write.table(t(Disease.input), file="$dir/Disease.input", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

chr=`cut -f 1 $dir/Disease.input`
region=`cut -f 1-3 $dir/Disease.input|sed 's/\t/:/;s/\t/-/;s/chr//'`

export PERL5LIB=software_needed/vcftools_0.1.11/lib/perl5/site_perl:$PERL5LIB
export PATH=software_needed/bin:$PATH
software_needed/vcftools_0.1.11/bin/vcf-query \
	-f "%CHROM:%POS\t%ID\n" Disease/00-All.vcf.gz $region > $dir/Disease.temp1
if [ ! -s $dir/Disease.temp1 ]; then exit; fi

for i in ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI
do
	myjoin -F 2 -f 2 $dir/Disease.temp1 Disease/ld_files/${chr}_${i}.txt |sed 's/$/\t'$i'/'| grep =  | cut -f 2,4,5,6,7 
done > $dir/Disease.temp2

if [ ! -s $dir/Disease.temp2 ]; then exit; fi

R --slave <<EOF
options(stringsAsFactors=FALSE)
library("plyr")
rsquare_cutoff <- 0.5
raw <- read.table("$dir/Disease.temp2")
raw[["V4"]] <- as.numeric( sprintf("%.2f",  raw[["V4"]] ) )
raw <- subset( raw, V4 > rsquare_cutoff )
get_max <- function(x){
		x[ which.max(x[,4]), ] 
	}
raw <- ddply(raw, .(V2, V5), get_max)
transform_table <- function(x){
		c( unique(x[,2]), unique(x[,3]), paste(x[,4], collapse=","), paste(x[,5], collapse=","))
	}
raw <- ddply(raw, .(V1, V2,V3), transform_table)
colnames(raw) <- c("tag SNP", "linked SNP", "linkage(rsquare)", "population")
write.table(raw, "$dir/Disease.temp3", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
EOF

myjoin $dir/Disease.temp3 <( cut -f 5,6 Disease/gwasCatalog_130621.hg19.UCSC ) | grep -P "^=" | cut -f 2,3,4,5,7 > $dir/Disease.output
}

run_disease $@
