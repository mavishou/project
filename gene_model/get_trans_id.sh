#!/bin/bash
# get all the trans ID
cat gencode.v19.annotation.gtf | grep -v "^##" | awk '{print $12}' | sed -e 's/"//g' -e 's/;//' | sort | uniq > all_trans.txt

# get all the lncRNAs ID
cat gencode.v19.long_noncoding_RNAs.gtf | grep -v "^##" | awk '{print $12}' | sed -e 's/"//g' -e 's/;//' | sort | uniq > all_lncrnas.txt
