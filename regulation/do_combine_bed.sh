# combine peaks
cd wgEncodeAwgTfbsUniform
for i in `ls -1 *Peak.gz`; do paste <(zcat $i | cut -f 1-3) <(zcat $i | cut -f 1 | sed 's/.*/'$i'/') ; done > All_peaks.bed
cd ../
mv wgEncodeAwgTfbsUniform/All_peaks.bed ./

# tidy the information about ChIP-seq peak files
echo "cell treatment antibody" | sed 's/ /\t/g' > peak_files.info
sed 's/; /\t/g' wgEncodeAwgTfbsUniform/files.txt | cut -f 1,7,8,9 | sed 's/cell=//;s/treatment=//;s/antibody=//' >> peak_files.info
