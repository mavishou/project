#!/bin/bash

function run_fastx_clipper {
	sample=$1
	cmd="fastx_clipper -Q33 -a TCGTATGCCGTCTTCTGCTTG -l 20 -v -z -i ${sample}.fastq -o ./remove_adapter/${sample}_remove_adapter.fastq.gz"
	echo $cmd
	$cmd
}

while read s
do
	now_date=`date +%y%m%d`
	(echo "Running fastx_clipper for $s"
	time run_fastx_clipper $s) 2>&1 | tee ./remove_adapter/logs/${now_date}-${s}_fastx_clipper.log
done