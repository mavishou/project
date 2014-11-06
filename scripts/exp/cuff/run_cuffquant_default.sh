#!/bin/bash

function run_cuffquant {
	sample=$1
	cmd="cuffquant -p 10 -o ../cuff/${sample}/ V4expCaculate.combined.gtf ${sample}.bam"
	echo "Running cuffquant for $sample..."
	echo $cmd
	eval $cmd
	echo "Done!"
}

while read s
do
	now_date=`date +%y%m%d`
	mkdir -p ../cuff/$s
	(echo $now_date
	time run_cuffquant $s) 2>&1 | tee ../cuff/logs/${now_date}-${s}_cuffquant.log
	echo "============================================="
done