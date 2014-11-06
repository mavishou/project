#!bin/bash

clip_dir=/lustre/user/houm/projects/AnnoLnc/CLIP
pipe_dir=$clip_dir/pipeclip
anno_dir=$clip_dir/annotation
logs_dir=$anno_dir/logs
mkdir -p $logs_dir

function run_annotation {
	sample=$1
	series=$2
	clip_type=$3
	indir=$pipe_dir/$series/$sample
	outdir=$anno_dir/$series
	mkdir -p $outdir
	
	if [ $clip_type = HITS ];then
		infile1=$indir/${sample}_CrossLinkingSites.substitution.bed
		infile2=$indir/${sample}_CrossLinkingSites.insertion.bed
		infile3=$indir/${sample}_CrossLinkingSites.deletion.bed
		outprefix1=${sample}_substitution
		outprefix2=${sample}_insertion
		outprefix3=${sample}_deletion
		cmd1="annotation -i $infile1 -o $outprefix1 -d $outdir"
		cmd2="annotation -i $infile2 -o $outprefix2 -d $outdir"
		cmd3="annotation -i $infile3 -o $outprefix3 -d $outdir"
		echo "Running annotation for $sample..."
		echo $cmd1
		eval $cmd1
		echo $cmd2
		eval $cmd2
		echo $cmd3
		eval $cmd3
	else
		infile1=$indir/${sample}_CrossLinkingSites.bed
		outprefix=$sample
		cmd="annotation -i $infile -o $outprefix -d $outdir"
		echo "Running annotation for $sample..."
		echo $cmd
		eval $cmd
	fi
	echo "Done!"
}

while read s e a t
do
	now_date=`date +%y%m%d`
	(echo $now_date
	time run_annotation $s $e $t) 2>&1 | tee $logs_dir/${now_date}-${s}_annotation.log
	echo "================================================="
done