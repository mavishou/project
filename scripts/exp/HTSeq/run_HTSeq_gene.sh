#!/bin/bash

function run_htseq_count {
	sample=$1
	cmd="htseq-count -f bam -r pos -s no ./${sample}.bam /lustre/user/houm/projects/AnnoLnc/V4expCaculate.combined.gtf > ../HTSeq/${sample}_gene.txt"
	echo $cmd
	htseq-count -f bam -r pos -s no ./${sample}.bam /lustre/user/houm/projects/AnnoLnc/V4expCaculate.combined.gtf > ../HTSeq/${sample}_gene.txt
}

while read s
do
	now_date=`date +%y%m%d`
	(echo "Running htseq-count for ${s}..."
	time run_htseq_count $s) 2>&1 | tee ../HTSeq/logs/${now_date}-${s}_htseq_count.log
done