#!/bin/bash

# 本来想写一个脚本，一边下载一边就能计算，下载完了的就计算，可是下载的速度好快，都快下载了一半了，我还没写完这个脚本，而且挺复杂的，算了吧。。。
# prepare the sample list


flag=1
while [ $flag -eq 1 ]
do
	ls *.bai | sed 's/.bam.bai$//' > now_bai
	myjoin now_bai already_run_cufflinks | grep ^+ | cut -f 2 > need_to_run_cufflinks

	if [ -s need_to_run_cufflinks ]
		then
		num=`wc -l need_to_run_cufflinks`
		if [ $num -lt 3 ]
			then
			cat need_to_run_cufflinks | bash run_cuffquant_default.sh
		else
			every_file_num=$(($num/3+1))
			first_end=$every_file_num
			seceond_start=$(($first_end+1))
			seconde_end=$(($every_file_num*2))
			third_start=$(($seconde_end+1))

			sed -n '1,${first_end}p' need_to_run_cufflinks > need_to_run_cufflinks_1
			sed -n '$seceond_start,${seconde_end}p' need_to_run_cufflinks > need_to_run_cufflinks_2
			sed -n '$third_start,${num}p' need_to_run_cufflinks > need_to_run_cufflinks_3

			cat need_to_run_cufflinks_1 | bash run_cuffquant_default.sh

	else
		sleep 5m
	fi
done