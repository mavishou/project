function run_mergeTrans {
dir=$1

bedtools merge -n -s -nms -i <(sort -k 1,1 -nk 2,2 $dir/blat.all.bed12) |sort -nrk 5,5| awk -F "\t" '{print $1"\t"$2"\t"$3"\tloci_"NR"\t"$4"\t"$6}' > $dir/merged.bed

}

run_mergeTrans $@
