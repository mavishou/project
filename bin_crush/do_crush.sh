dir=$1

##. Cpding potential prediction
#if [ ! -e $dir/input.fa ] ; then echo "There is no $dir/input.fa"; exit; fi
#bash bin_crush/do_CPC.sh $dir &
#
##. Predict RNA secondary structure
#if [ ! -e $dir/input.fa ] ; then echo "There is no $dir/input.fa"; exit; fi
#bash bin_crush/do_SecStructure.sh $dir &
#
##. Map input sequence onto hg19 using blat
#if [ ! -e $dir/input.fa ] ; then echo "There is no $dir/input.fa"; exit; fi
#bash bin_crush/do_blat.sh $dir 

##. Check the number of locis, only one loci at a time
#if [ ! -e $dir/ExpCall.gtf ]; then echo "There is no $dir/ExpCall.gtf"; exit ; fi
#loci_check=`bash bin_crush/single_loci_check.sh $dir`
#if [ $loci_check -ne 0 ]; then echo "didn't pass single loci check "; exit; fi

#. Combine with known genes
if [ ! -e $dir/ExpCall.gtf ]; then echo "There is no $dir/ExpCall.gtf"; exit ; fi
if [ ! -e $dir/exons.bed ]; then echo "There is no $dir/exons.bed"; exit ; fi
bash bin_crush/do_combineKnown.sh $dir

#. Characterizing expression regulation context
if [ ! -e $dir/exons.bed ]; then echo "There is no $dir/exons.bed"; exit; fi
bash bin_crush/do_TransReg.sh $dir 3000 1000 &

#. Calculate evolution score
if [ ! -e $dir/exons.bed ]; then echo "There is no $dir/exons.bed"; exit; fi
bash bin_crush/do_Evolution.sh $dir &

##. Disease
#if [ ! -e $dir/exons.bed ]; then echo "There is no $dir/exons.bed"; exit; fi
#bash bin_crush/do_Disease.sh $dir 3000 1000 &

. Calculate FPKM using cufflinks
if [ ! -e $dir/exons.bed ]; then echo "There is no $dir/exons.bed"; exit; fi
bash bin_crush/do_ExpCall.sh $dir

#. GO annotation 
if [ ! -e $dir/ExpCall.output ]; then echo "There is no $dir/ExpCall.output"; exit ; fi
bash bin_crush/do_GOanno.sh $dir
