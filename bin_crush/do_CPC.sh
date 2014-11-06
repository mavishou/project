function run_CPC {
dir=$1

if [ -d $dir/CPC.temp ]; then rm -fr $dir/CPC.temp; fi
mkdir $dir/CPC.temp
if [ -e $dir/CPC.output ]; then rm -fr $dir/CPC.output; fi

bash software_needed/CPC/bin/run_predict_local.sh $dir/input.fa $dir/CPC.output $dir/CPC.temp/ $dir/CPC.temp/my

rm -rf $dir/CPC.temp
}

run_CPC $@
