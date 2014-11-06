#!/bin/bash

function run_cuffquant {
	sample=$1
	cmd="cuffquant -p 20 --no-effective-length-correction -o ../cuff/effective2/${sample}/  -u /lustre/user/houm/projects/AnnoLnc/V4expCaculate.combined.gtf ${sample}.bam"
	echo $cmd
	eval $cmd
}

while read s
do
	now_date=`date +%y%m%d`
	mkdir -p ../cuff/effective/$s
	(echo $now_date
	echo "Running cuffquant for $s"
	time run_cuffquant $s) 2>&1 | tee ../cuff/effective2/logs/${now_date}-${s}_cuffquant_effective.log
done