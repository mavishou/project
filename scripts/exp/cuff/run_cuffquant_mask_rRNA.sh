#!/bin/bash

function run_cuffquant {
	sample=$1
	cmd="cuffquant -p 10 -o ../cuff/mask_rRNA/${sample}/ -M /lustre/user/houm/projects/AnnoLnc/expression/rRNA.gtf  /lustre/user/houm/projects/AnnoLnc/V4expCaculate.combined.gtf ${sample}.bam"
	echo $cmd
	$cmd
}

while read s
do
	now_date=`date +%y%m%d`
	mkdir -p ../cuff/mask_rRNA/$s
	(echo "Running cuffquant for $s"
	time run_cuffquant $s) 2>&1 | tee ../cuff/mask_rRNA/logs/${now_date}-${s}_cuffquant_mask_rRNA.log
done