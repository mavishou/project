#!/bin/bash
# stat the stat_file

# the stat_file is from JangS
(stat_file=gencode_whole_lncRNAdb.final.format
total_records=`wc -l $stat_file | awk '{print $1}'`
mapped_gencode=`cut -f 1 $stat_file | sort | uniq | wc -l`
mapped_lncrnadb=`cut -f 2 $stat_file | sort | uniq | wc -l`

echo "# Total records: $total_records"
echo "# Mapped GENCODE lncRNAs: $mapped_gencode"
echo "# Mapped lncrnadb lncRNAs: $mapped_lncrnadb"
echo
echo "---------------------------------------------"
echo

echo "Multi-mapped GENCODE lncRNAs:"
cut -f 1 gencode_whole_lncRNAdb.final.format | sort | uniq -c | awk '$1>1 {print $1"\t"$2}'
echo
echo "Multi-mapped lncrnadb lncRNAs:"
cut -f 2 gencode_whole_lncRNAdb.final.format | sort | uniq -c | awk '$1>1 {print $1"\t"$2}'
echo 
echo "---------------------------------------------"
echo) > gencode_lncrnadb_map.stat
