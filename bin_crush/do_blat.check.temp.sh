function run_blatCheck {
dir=$1
if [ $dir = "" ]; then echo "Please input Working Directory"; exit; fi

if [ -e $dir/test.fa ] ; then rm -f $dir/test.fa ; fi
if [ -e $dir/input.fa ] ; then rm -f $dir/input.fa ; fi

gtf_to_fasta $dir/test.gtf blat/hg19.fa $dir/test.fa
strand=`cut -f 7 $dir/test.gtf | sort -u `
if [ $strand = "-" ]; then revseq $dir/test.fa $dir/input.fa ;fi
if [ $strand = "+" ]; then mv $dir/test.fa $dir/input.fa ;fi

bash bin_crush/do_blat2.sh $dir 

}

run_blatCheck $@
