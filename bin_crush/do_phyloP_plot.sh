function run_phyloP_plot {
dir=$1
loci=$2

for id in `myjoin <(grep $loci $dir/merged.bed | cut -f 5 | sed 's/,/\n/g') $dir/SecStructure.output/IDs| grep = | cut -f 2`
do
grep -w "$id" $dir/blat.psl > $dir/pos_map.psl
if [ -e $dir/pos_map.psl ] && [ -e ${dir}/SecStructure.output/${id}.seq.txt ] && [ -e ${dir}/SecStructure.output/${id}.struct.txt ] ; then 
seq=`cat ${dir}/SecStructure.output/${id}.seq.txt`
struct=`cat ${dir}/SecStructure.output/${id}.struct.txt`
for p in Primates Placental All
do
R --slave <<EOF
options(stringsAsFactors=FALSE)
scores <- scan("${dir}/Evolution.output/phyloP.${id}.${p}.nt_score", what="numeric")
pos_map <- read.table("$dir/pos_map.psl", sep="\t")
trans_length <- as.numeric( pos_map[1,11] )
strand <- pos_map[1,9]
block_size <- as.numeric( unlist( strsplit(pos_map[1,19], ",") ) )
Q_start <- as.numeric( unlist( strsplit( pos_map[1,20], ",") ) )
T_start <- as.numeric( unlist( strsplit( pos_map[1,21], ",") ) )
trans_pos <- vector()
for( i in 1:length(block_size)){
trans_pos <- c( trans_pos, seq( from=Q_start[i], by=1, length.out=block_size[i] ) )
}
new_scores <- rep(0, trans_length)
new_scores[trans_pos+1] <- scores
if(strand %in% c("-", "0")){ new_scores <- rev(new_scores) }
write( new_scores, file="${dir}/Evolution.output/phyloP.${id}.${p}.new.nt_score", ncolumns=1)
EOF
phyloP=`cat ${dir}/Evolution.output/phyloP.${id}.${p}.new.nt_score | tr "\n" ";"`
java -Djava.awt.headless=true -cp software_needed/bin/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd -colorMapMax 2  -colorMapMin -2 -title "secondary structure of ${id} \n color by ${p} phyloP Score" -algorithm naview -sequenceDBN $seq -structureDBN $struct -colorMap $phyloP -o ${dir}/SecStructure.output/phyloP.${id}.${p}.EPS
convert -density 150 -geometry 100% ${dir}/SecStructure.output/phyloP.${id}.${p}.EPS ${dir}/SecStructure.output/phyloP.${id}.${p}.png
done
fi
done
}
run_phyloP_plot $@
