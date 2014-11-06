#!/bin/bash
(zcat remove_adaptor_3.fa.gz | fastx_collapser | gzip > collapsed_3.fa.gz
bowtie /lustre/user/houm/genome/human/hg19/hg19 -a -p 10 -v 2 -m 10 --best --strata -t -f <(zcat collapsed_3.fa.gz) > collapsed_3.bt) 2>&1 | tee run_140306.log
