#!bin/bash

clip_dir=/lustre/user/houm/projects/AnnoLnc/CLIP
bwa_dir=$clip_dir/bwa
pipe_dir=$clip_dir/pipeclip
logs_dir=$pipe_dir/logs
mkdir -p $logs_dir

function run_pipeclip {
	sample=$1
	series=$2
	clip_type=$3
	indir=${bwa_dir}/$series
	outdir=$pipe_dir/$series
	mkdir -p $outdir
	infile=$indir/${sample}_unique.bam
	
	clip_code=1
	if [ $clip_type = HITS ];then
		clip_code=0
	fi
	if [ $clip_type = iCLIP ];then
		clip_code=3
	fi
	if [ $clip_type = PAR6 ];then
		clip_code=2
	fi

	cmd="pipeclip -i $infile -c $clip_code -f 0.05 -F 0.05 -o $sample -d $outdir"
	echo "Running piepclip for $sample..."
	echo $cmd
	eval $cmd
	echo "Done!"
}

while read s e a t
do
	now_date=`date +%y%m%d`
	(echo $now_date
	time run_pipeclip $s $e $t) 2>&1 | tee $logs_dir/${now_date}-${s}_pipeclip.log
	echo "================================================="
done