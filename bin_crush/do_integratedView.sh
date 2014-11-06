function run_integratedView {
dir=$1
now_dir=`pwd`

if [ -e $dir/Disease.output ]; then  grep -wf <(cut -f 2 $dir/Disease.output| sort -u ) $dir/Disease.temp1 | sed 's/^.*://' | sort -u > $dir/Disease.plot ; fi
#input files: Evolution.input TransReg.temp Evolution.output/ Disease.plot
R --slave <<EOF
#library("Gviz")
#library("GenomicRanges")
library(RSclient)
conn<-RS.connect(host="127.0.0.1", port=6312)
options(stringsAsFactors=FALSE)

tryCatch.my <- function(expr){
tryCatch(expr, error = function(e)e )
}

region <- read.table("$dir/Evolution.input", sep="\t")
from <- region[,2]
to <-region[,3]
chr <- region[,1]
RS.assign(conn, "to", to)
RS.assign(conn, "from", from)
RS.assign(conn, "chr", chr)

bed <- read.table("$dir/exons.bed", sep="\t")
RS.assign(conn, "bed", bed)
trans_id <- unique(bed[["V4"]])
RS.assign(conn, "trans_id", trans_id)
grTrack <- RS.eval(conn, GeneRegionTrack( start=min(bed[["V2"]]), end=max(bed[["V3"]]), rstart=bed[["V2"]], rends=bed[["V3"]], chromosome=chr, genome="hg19", transcript=bed[["V4"]], symbol=bed[["V4"]], name="Gene\nModel", strand=bed[["V6"]], feature="exon", showId=TRUE, fontcolor="black" ) )

files <- read.table("TransReg/peak_files.info")
e <- tryCatch.my( overlaps <- read.table("$dir/TransReg.temp") )
if( is.null( e[["message"]] ) ) { 
ChIP_exp <- overlaps[, 10 ]
overlaps[, 10 ] <- paste( files[ ChIP_exp, 1], files[ChIP_exp,3], sep=', ')
overlaps <- unique( overlaps[, c("V8","V9","V10") ] )
overlaps <- overlaps[ 1:min(100, nrow(overlaps)), ]
gr <- overlaps[, "V10" ] 
RS.assign(conn, "gr", gr)
RS.assign(conn, "overlaps", overlaps)
TFBS <- RS.eval(conn, AnnotationTrack(start=overlaps[,"V8"], end=overlaps[,"V9"], strand="*", chromosome=chr, genome="hg19", feature="TFBS",id=gr, name="TFBS", stacking="full", fontcolor="black", fontsize=3) ) } 

conservation <- list()
for(i in c("DAF.wig", "phyloP.Primates.wig", "phyloP.Placental.wig", "phyloP.All.wig")){
	e <- tryCatch.my( conservation[[i]] <- read.table(paste("$dir/Evolution.output/", i, sep=""), sep="\t" ) )
}

RS.assign(conn, "conservation", conservation)
if( ! is.null( conservation[["DAF.wig"]] ) ){
DAF <- RS.eval(conn, DataTrack(genome = "hg19", chromosome = chr, width=1, start=conservation[["DAF.wig"]][,1], data=conservation[["DAF.wig"]][,2], type = "hist", col.histogram = "darkblue", fill.histogram = "darkblue", ylim = c(0, 1), name = "DAF") ) } 

if( ! is.null( conservation[["phyloP.Primates.wig"]] ) ){
phyloP_primates <- RS.eval(conn, DataTrack(genome = "hg19", chromosome = chr, width=1, start=conservation[["phyloP.Primates.wig"]][,1], data=conservation[["phyloP.Primates.wig"]][,2], type = "hist", col.histogram = "darkblue", fill.histogram = "darkblue", ylim = c(-3.7, 4), name = "phyloP\nprimates", window = "auto") ) } 
if( ! is.null( conservation[["phyloP.Placental.wig"]] ) ){
phyloP_placental <- RS.eval(conn, DataTrack(genome = "hg19", chromosome = chr, width=1, start=conservation[["phyloP.Placental.wig"]][,1], data=conservation[["phyloP.Placental.wig"]][,2], type = "hist", col.histogram = "darkblue", fill.histogram = "darkblue", ylim = c(-3.7, 4), name = "phyloP\nmammals", window = "auto") ) } 
if( ! is.null( conservation[["phyloP.All.wig"]] ) ){
phyloP_vertebrates <- RS.eval(conn, DataTrack(genome = "hg19", chromosome = chr, width=1, start=conservation[["phyloP.All.wig"]][,1], data=conservation[["phyloP.All.wig"]][,2], type = "hist", ylim = c(-3.7, 4), col.histogram = "darkblue", fill.histogram = "darkblue", name = "phyloP\nvertebrates", window = "auto") ) } 

e <- tryCatch.my( SNP <- subset( read.table("$dir/Disease.plot", sep="\t"), V1 <= to & V1 >= from ) )
if( is.null( e[["message"]] ) & "SNP" %in% ls() ){
if( nrow(SNP) > 0 ){
DiseaseTracks <- list()
for( i in 1:nrow(SNP) ){
RS.assign(conn, "SNP", SNP[i,])
DiseaseTracks[[ SNP[i,2] ]] <- RS.eval(conn, AnnotationTrack(genome = "hg19", chromosome = chr, width=1, start=SNP[,1], id=SNP[,2], fontcolor="black", fill = "red", name = "DiseaseSNP", fontsize=8, shape = "box") ) } }
}

axTrack <- RS.eval(conn, GenomeAxisTrack())
idxTrack <- RS.eval(conn, IdeogramTrack(genome = "hg19", chromosome = chr) )

plot_tracks <- c("idxTrack","axTrack","grTrack","TFBS","DAF","phyloP_primates","phyloP_placental","phyloP_vertebrates")
plot_tracks <- plot_tracks[ plot_tracks %in% ls() ]

RS.eval(conn, jpeg("$now_dir/$dir/GB_plot.jpeg", res=200, height=1500, width=2500) )
plot_tracks <- mget(plot_tracks)
if( "DiseaseTracks" %in% ls()){ plot_tracks <- c(plot_tracks, DiseaseTracks) }
RS.assign(conn, "plot_tracks", plot_tracks)
RS.eval(conn, plotTracks( plot_tracks, from = from, to = to, showTitle = TRUE, showFeatureId=TRUE, background.title = "darkblue") )
RS.eval(conn, dev.off() )
RS.close(conn)
EOF
#rm -f $dir/Disease.plot
}

run_integratedView $@
