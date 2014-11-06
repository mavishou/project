function run_blat {
dir=$1

if [ -d $dir/inDir ]; then rm -rf $dir/inDir/; fi
mkdir $dir/inDir/
software_needed/blat/gfClient -t=dna -q=rna 127.0.0.1 33333 ./ $dir/input.fa $dir/inDir/blat.psl

if [ -d $dir/temp ]; then rm -rf $dir/temp/; fi
if [ -e $dir/sorted.psl ]; then rm -f $dir/sorted.psl; fi
pslSort dirs $dir/sorted.psl $dir/temp $dir/inDir
if [ -e $dir/blat.psl ]; then rm -f $dir/blat.psl; fi
pslReps $dir/sorted.psl  $dir/blat.psl $dir/out.psr

if [ -e $dir/blat.all.bed12 ]; then rm -f $dir/blat.all.bed12; fi
perl bin_crush/psl2bed.pl $dir/blat.psl $dir/blat.all.bed12

if [ ! -s $dir/blat.all.bed12 ]; then exit 1 ; fi

# This bed is 1 based
if [ -e $dir/exons.all.bed ]; then rm -f $dir/exons.all.bed; fi
bedtools bed12tobed6 -i $dir/blat.all.bed12 | awk -F "\t" '{print $1"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6}' > $dir/exons.all.bed
#bed12tobed6 is wrong when there  is only one exon
R --slave <<EOF
options(stringsAsFactors=FALSE)
bed <- read.table("$dir/exons.all.bed", sep="\t")
trans <- bed[,4]
singleExonTrans <- levels( as.factor(trans) )[ table(trans) == 1 ] 
bed_new <- bed
for( i in 1:nrow(bed_new) ){
	if( bed_new[ i, 4] %in% singleExonTrans) { bed_new[ i, 3] <- bed_new[ i, 3] - 1 }
}
write.table(bed_new, file="$dir/exons.all.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

if [ -e $dir/ExpCall.all.gtf ]; then rm -f $dir/ExpCall.all.gtf; fi
bash bin_crush/bed2gtf.sh $dir/exons.all.bed $dir/ExpCall.all.gtf 

}

run_blat $@
