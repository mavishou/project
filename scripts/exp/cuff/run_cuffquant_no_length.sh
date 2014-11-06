#!/bin/bash

function run_cuffquant {
	sample=$1
	cmd="cuffquant -p 10 -o ../cuff/no_length/${sample}/ -M /lustre/user/houm/projects/AnnoLnc/expression/rRNA.gtf --no-length-correction /lustre/user/houm/projects/AnnoLnc/V4expCaculate.combined.gtf ${sample}.bam"
	echo $cmd
	$cmd
}

while read s
do
	now_date=`date +%y%m%d`
	mkdir -p ../cuff/no_length/$s
	(echo "Running cuffquant for $s"
	time run_cuffquant $s) 2>&1 | tee ../cuff/no_length/logs/${now_date}-${s}_cuffquant_no_length.log
done