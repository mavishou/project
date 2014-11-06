#!/bin/bash
(zcat SRR764667.fastq.gz SRR764668.fastq.gz SRR764669.fastq.gz | fastq_to_fasta -r -v -Q 33 |  fastx_clipper -a TCGTATGCCGTCTTCTGCTTG -l 20 -v |  fastx_collapser | gzip > LIN28B_total.fa.gz
bowtie /lustre/user/houm/genome/human/hg19/hg19 -a -p 10 -v 2 -m 10 --best --strata -t -f <(zcat LIN28B_total.fa.gz) > LIN28B_total.bt 
) 2>&1 | tee run_LIN28B_140307.log
