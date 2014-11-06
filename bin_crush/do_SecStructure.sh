function run_SecStructure {
dir=$1

now_dir=`pwd`
cd $dir
if [ -d SecStructure.output ]; then rm -rf SecStructure.output; fi
mkdir SecStructure.output
cd SecStructure.output

cat ../input.fa| perl -e '$t_num=1;while(<STDIN>){$line=$_; if($line =~ /^>/){print "\n>$t_num\n"; $t_num++; } else {chomp($line); print $line;} }'| $now_dir/software_needed/ViennaRNA/bin/RNAfold -p -d2 -noLP | grep -P -A 2 "^>" | sed 's/(-[0-9]*\.[0-9]*)//' > seq_struct.txt 
grep -P "^>" ../input.fa | sed 's/>//' > IDs 

for INP in `ls *_ss.ps`
do
	id=`basename $INP _ss.ps`
	newname=`head -n $id IDs|tail -n 1`
	#$now_dir/software_needed/ViennaRNA/ViennaRNA/bin/relplot.pl $INP ${id}_dp.ps > ${id}_rss.ps
	#convert -density 300 -geometry 100% ${id}_rss.ps $newname.png
	grep -P -A 2 "^>${id}" seq_struct.txt | head -n 2 | tail -n 1 > ${newname}.seq.txt 
	grep -P -A 2 "^>${id}" seq_struct.txt | head -n 3 | tail -n 1 > ${newname}.struct.txt 
	$now_dir/software_needed/ViennaRNA/ViennaRNA/bin/relplot.revised.pl $INP ${id}_dp.ps > ${newname}.entropy.txt
	seq=`cat ${newname}.seq.txt`
	struct=`cat ${newname}.struct.txt`
	entropy=`cat ${newname}.entropy.txt`
java -Djava.awt.headless=true -cp $now_dir/software_needed/bin/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd  -colorMapMax 1  -colorMapMin 0 -title "secondary structure of ${newname} \n color by entropy (RNAfold)" -algorithm naview -sequenceDBN $seq -structureDBN $struct -colorMap $entropy -o ${newname}.entropy.EPS
convert -density 150 -geometry 100% ${newname}.entropy.EPS ${newname}.entropy.png
done
#rm -f *.ps
cd $now_dir
}

run_SecStructure $@
