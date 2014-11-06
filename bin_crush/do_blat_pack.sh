function run_blat_pack {
dir=$1

#. Map input sequence onto hg19 using blat
if [ ! -e $dir/input.fa ] ; then echo "There is no $dir/input.fa"; exit; fi
bash bin_crush/do_blat2.sh $dir
if [ ! -e $dir/exons.all.bed ] ; then echo "There is no $dir/exons.all.bed"; exit; fi
if [ ! -e $dir/ExpCall.all.gtf ] ; then echo "There is no $dir/ExpCall.all.gtf"; exit; fi
bash bin_crush/do_combineKnown.sh $dir 

if [ ! -e $dir/merged.bed ] ; then echo "There is no $dir/merged.bed"; exit; fi
while read c s e loci trans s 
do
#Prepare blat.bed12, exons.bed and ExpCall.gtf
if [ -e $dir/blat.bed12 ] ; then rm -f $dir/blat.bed12 ; fi
bedtools intersect -a $dir/blat.all.bed12 -b <(grep $loci $dir/merged.bed) -s -wa  -u -split > $dir/blat.bed12
if [ -e $dir/exons.bed ] ; then rm -f $dir/exons.bed ; fi
bedtools intersect -a $dir/exons.all.bed -b <(grep $loci $dir/merged.bed) -s -wa  -u > $dir/exons.bed
if [ -e $dir/ExpCall.gtf ] ; then rm -f $dir/ExpCall.gtf ; fi
bedtools intersect -a $dir/ExpCall.all.gtf -b <(grep $loci $dir/merged.bed) -s -wa  -u > $dir/ExpCall.gtf

if [ -e $dir/miRBS.${loci}.out ] ; then rm -rf $dir/miRBS.${loci}.out ; fi
if [ -e $dir/TransReg.${loci}.out ] ; then rm -rf $dir/TransReg.${loci}.out ; fi
if [ -e $dir/Disease.${loci}.out ] ; then rm -rf $dir/Disease.${loci}.out ; fi
if [ -d $dir/Evolution.${loci}.out ] ; then rm -rf $dir/Evolution.${loci}.out ; fi
if [ -d $dir/GB_plot.${loci}.jpeg ] ; then rm $dir/Evolution.${loci}.jpeg ; fi

#miRNA binding sites prediction
if [ ! -e $dir/blat.bed12 ] ; then echo "There is no $dir/blat.bed12" ;  exit; fi
bash bin_crush/do_miRBS.sh $dir
#transcriptional regulation
if [ ! -e $dir/exons.bed ] ; then echo "There is no $dir/exons.bed" ;  exit; fi
bash bin_crush/do_TransReg.sh $dir 5000 1000
#Disease
if [ ! -e $dir/exons.bed ] ; then echo "There is no $dir/exons.bed" ;  exit; fi
bash bin_crush/do_Disease.sh $dir 5000 1000
#Evolution
if [ ! -e $dir/exons.bed ] ; then echo "There is no $dir/exons.bed" ;  exit; fi
bash bin_crush/do_Evolution.hgwiggle.sh $dir 5000 1000
#do_integratedView.sh input conservation scores from do_Evolution.hgwiggle.sh
if [ ! -e $dir/exons.bed ] ; then echo "There is no $dir/exons.bed" ;  exit; fi
bash bin_crush/do_integratedView.sh $dir
#color RNA secondary structure with phyloP score
if [ ! -e $dir/blat.psl ] ; then echo "There is no $dir/blat.psl" ;  exit; fi
bash bin_crush/do_phyloP_plot.sh $dir $loci


if [ -e $dir/miRBS.output ] ; then mv $dir/miRBS.output $dir/miRBS.${loci}.out ; fi
if [ -e $dir/TransReg.output ] ; then mv $dir/TransReg.output $dir/TransReg.${loci}.out ; fi
if [ -e $dir/Disease.output ] ; then mv $dir/Disease.output $dir/Disease.${loci}.out ; fi
if [ -d $dir/Evolution.output ] ; then mv $dir/Evolution.output $dir/Evolution.${loci}.out ; fi
if [ -e $dir/GB_plot.jpeg ] ; then mv $dir/GB_plot.jpeg $dir/GB_plot.${loci}.jpeg ; fi

done < $dir/merged.bed
}

run_blat_pack $@
