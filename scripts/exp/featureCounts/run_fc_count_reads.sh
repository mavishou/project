#!/bin/bash

now_date=`date +%y%m%d`

sample_list_file=$1
list=`sed 's/$/.bam/' $sample_list_file | tr '\n' ' '`
cmd="featureCounts -T 20 -t exon -g gene_id -a V4expCaculate.combined.gtf -o ../featureCounts/32_sample_gene_reads_count.txt $list"
(echo $now_date
echo "Running featureCounts......"
echo $cmd
time eval $cmd) 2>&1 | tee ../featureCounts/logs/${now_date}-32_samples-featureCouts.log
