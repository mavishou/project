#!/bin/bash

outdir=/lustre/user/houm/projects/AnnoLnc/CLIP/bwa
genome=/lustre/user/houm/genome/human/hg19/hg19.fa
header=/lustre/user/houm/genome/human/hg19/hg19_sam_header

function run_bwa {
	sample=$1
	sai_file=${outdir}/${sample}.sai
	fastq_file=${sample}_remove_adapter.fastq
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
	echo "Done"

	if [ -e $bam_file ];then
		rm -f $sam_file
	fi

	if [ -e $unique_bam_file ];then
		rm -f $unique_sam_file
	fi
}

while read s
do
	now_date=`date +%y%m%d`
	(echo $now_date
	time run_bwa $s
	echo Done!
	echo -e "======================================================\n") 2>&1 | tee ${outdir}/logs/${now_date}-${s}_bwa.log
done