function run_single_loci_check {
dir=$1
if [ ! -e $dir/ExpCall.gtf ]; then echo "There is ni $dir/ExpCall.gtf"; exit; fi

export PATH=software_needed/cufflinks-2.0.1.Linux_x86_64:$PATH
loci_check=`cuffmerge -o $dir/merged_asm -g $dir/ExpCall.gtf <( echo $dir/ExpCall.gtf ) 2>&1 | grep "Processed 1 loci" | wc -l `
rm -rf $dir/merged_asm
if [ $loci_check -eq 0 ]; then echo 1; exit; fi
echo 0
}

run_single_loci_check $@
