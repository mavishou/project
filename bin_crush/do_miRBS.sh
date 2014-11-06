function run_miRBS {
dir=$1

if [ ! -e $dir/blat.bed12 ] ; then echo "There is no $dir/blat.bed12" ;  exit 1; fi
if [ -e $dir/ma.fa ] ; then rm -f $dir/ma.fa ; fi
if [ -e $dir/targetscan.out ] ; then rm -f $dir/targetscan.out ; fi
if [ -e $dir/miRBS.out ] ; then rm -f $dir/miRBS.out ; fi

#Extract alignment from UCSC 46-way alignment
chr=`cut -f 1 ${dir}/blat.bed12| sort -u`
export PYTHONPATH=/home/tangx/software/galaxy-galaxy-dist-5c789ab4144a/lib/
/home/tangx/galaxy-python/python /home/tangx/software/galaxy-galaxy-dist-5c789ab4144a/tools/maf/interval_maf_to_merged_fasta.py -d hg19 -i $dir/blat.bed12 -m miRBS/maf/${chr}.maf -o $dir/ma.fa -G -I miRBS/maf/${chr}.maf.index --overwrite_with_gaps=TRUE --mafSourceType=user

#Predict miRNAs' binding sites by targetScan
bin_crush/myjoin miRBS/species2IDs <(cat $dir/ma.fa | paste - - | sed 's/^>//' | sed 's/\./\t/' ) | grep "^=" | awk -F "\t" '{print $5"\t"$3"\t"$6}' > $dir/ma_seq.txt
perl software_needed/targetscan/targetscan_60.pl miRBS/miR_for_targetscan_human.miRCode.txt $dir/ma_seq.txt $dir/targetscan.out
awk -F "\t" '{if($3==9606)print}' $dir/targetscan.out | cut -f 1,2,6,7,9,12 > $dir/miRBS.out

#Convert species' ID to names
R --slave <<EOF
options(stringsAsFactors=FALSE)
library("plyr")

tryCatch.my <- function(expr){
tryCatch(expr, error = function(e)e )
}

e <- tryCatch.my(raw <- read.table(file="$dir/miRBS.out", sep="\t") )
if( is.null( e[["message"]] ) ) {

species2IDs <- read.table(file="miRBS/species2IDs", sep="\t", row.names=2)
species2lineages <- read.table(file="miRBS/species2lineages", sep="\t", row.names=1)
t <- as.factor( as.vector( as.matrix( species2lineages) ) )
names(t) <- rownames(species2lineages)
species2lineages <- t

convert_id2name <- function(x){ paste( species2IDs[ as.character( x), 1 ], collapse=",") }
group_lineage <- function(x){ x <- unlist( strsplit(x, " ") ); 
			tapply( x, species2lineages[as.character(x)], convert_id2name )	}
lineage_score <- function(x){ x <- unlist( strsplit(x, " ") ); 
			tapply( x, species2lineages[as.character(x)], length ) / table(species2lineages) }

groups <- laply( raw[,6], group_lineage, .drop = FALSE)[,c("Primates", "Mammals", "Vertebrates"), drop = FALSE]
scores <- laply( raw[,6], lineage_score, .drop = FALSE)[,c("Primates", "Mammals", "Vertebrates"), drop = FALSE]
scores <- apply( scores, c(1,2), function(x){ if(is.na(x)){ x<-0} ; sprintf("%.2f", x) } )
new <- cbind( raw[, 1:5], scores )
new <- cbind( new, groups )
colnames(new) <- c("lncRNA", "miRNA family", "match start", "match end", "match type", "conservation score (primates)", "conservation score (non-primate mammals)", "conservatoin score (non-mammal vertebrates)", "matched species (primates)", "matched species (non-primate mammals)", "matched species (non-mammal vertebrates)" )

write.table(new, file="$dir/miRBS.output", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

}
EOF
}

run_miRBS $@
