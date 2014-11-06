#!/bin/bash

clip_dir=/lustre/user/houm/projects/AnnoLnc/CLIP
fq_dir=$clip_dir/fastq
bwa_dir=$clip_dir/bwa
logs_dir=$bwa_dir/logs
mkdir -p $logs_dir
genome=/lustre/user/houm/genome/human/hg19/hg19.fa
header=/lustre/user/houm/genome/human/hg19/hg19_sam_header

function run_bwa {
	sample=$1
	series=$2
	outdir=$bwa_dir/$series
	mkdir -p $outdir
	indir=$fq_dir/$series
	sai_file=${outdir}/${sample}.sai
	fastq_file=$indir/${sample}_remove_adapter.fastq
	sam_file=${outdir}/${sample}.sam
	bam_file=${outdir}/${sample}.bam
	unique_sam_file=${outdir}/${sample}_unique.sam
	unique_bam_file=${outdir}/${sample}_unique.bam

	cmd1="bwa aln -n 1 -i 0 -t 20 $genome $fastq_file > $sai_file"
	echo "Running bwa aln for ${sample}..."
	echo $cmd1
	eval $cmd1
	echo "Done!"
	echo
	cmd2="bwa samse $genome $sai_file $fastq_file > $sam_file"
	echo "Running bwa samse for ${sample}..."
	echo $cmd2
	eval $cmd2
	echo "Done!"	
	echo
	echo "Filtering unique mapped reads..."
	cmd3="samtools view -S -F 4 $sam_file | grep -w 'XT:A:U' > $unique_sam_file"
	echo $cmd3
	eval $cmd3
	echo "Done!"
	echo
	echo -n "Converting sam to bam..."
	cmd4="samtools view -Sbh <(cat $header $unique_sam_file) > $unique_bam_file"
	echo $cmd4
	eval $cmd4
	cmd5="samtools view -Sbh $sam_file > $bam_file"
	echo $cmd5
	eval $cmd5
	echo "Done!"

	if [ -e $bam_file ];then
		rm -f $sam_file
	fi

	if [ -e $unique_bam_file ];then
		rm -f $unique_sam_file
	fi
}

while read s e
do
	now_date=`date +%y%m%d`
	(echo $now_date
	time run_bwa $s $e
	echo -e "======================================================\n") 2>&1 | tee $logs_dir/${now_date}-${s}_bwa.log
done