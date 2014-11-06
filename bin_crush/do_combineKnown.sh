function run_combineKnown {
dir=$1

export PATH=software_needed/cufflinks-2.0.1.Linux_x86_64:$PATH
if [ -d $dir/cuffcompare ]; then rm -rf $dir/cuffcompare; fi
mkdir $dir/cuffcompare
new_dir=`pwd`
cd $dir/cuffcompare
cuffcompare -T -r ../ExpCall.all.gtf $new_dir/combineKnown/V4expCaculate.combined.gtf 
awk -F "\t" '{if($4=="j") print $5}' cuffcmp.tracking | sed 's/|/\t/g' | cut -f 2 | sort -u > overlapped.ref.IDs

cuffcompare -T -Rr $new_dir/combineKnown/V4expCaculate.combined.gtf ../ExpCall.all.gtf
paste <(cut -f 2 cuffcmp.loci | sed 's/\[/\t/; s/\]/\t/' |cut -f 1,3,4 | sed 's/-/\t/') <( cut -f 1 cuffcmp.loci) <( awk -F "\t" '{print $4","$3}' cuffcmp.loci| sed 's/|/=/g' | sed 's/,-$//') <(cut -f 2 cuffcmp.loci | sed 's/\[/\t/; s/\]/\t/' |cut -f 2) > merged.bed
R --slave <<EOF
options(stringsAsFactors=FALSE)
library('plyr')
bed <- read.table("merged.bed", sep="\t")
extract_trans <- function(x){ x <- unlist( strsplit(x, ",") ) ; paste( sub(".*=", "", x), collapse="," )}
bed[,5] <- unlist( llply( bed[,5], extract_trans ) ) 
write.table( bed , "merged.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
EOF
mv merged.bed ../

cd ../
if [ ! -s cuffcompare/overlapped.ref.IDs ]; then exit; fi
grep -f cuffcompare/overlapped.ref.IDs $new_dir/combineKnown/V4expCaculate.combined.gtf >> ExpCall.all.gtf
#rm -rf cuffcompare
sed 's/; /\t/g' ExpCall.all.gtf| awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$10"\t.\t"$7}' | sed 's/transcript_id //' > exons.all.bed
cd $new_dir
}
run_combineKnown $@
