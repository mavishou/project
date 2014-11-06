cd /rd1/user/liufl/software/apache-tomcat-6.0.37/working_crush/C2797E20-76BF-11E3-958C-E4B401005D36/SecStructure.output
newname=lincRNA-RoR
seq=`cat ${newname}.seq.txt`
struct=`cat ${newname}.struct.txt`
entropy=`cat ${newname}.entropy.txt`
now_dir=/rd1/user/tangx/crack

java -Djava.awt.headless=true -cp $now_dir/software_needed/bin/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd  -colorMapMax 1  -colorMapMin 0 -title "secondary structure of ${newname} \n color by entropy (RNAfold)" -algorithm naview -sequenceDBN $seq -structureDBN $struct -colorMap $entropy -o ${newname}.entropy.svg

java -Djava.awt.headless=true -cp $now_dir/software_needed/bin/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd  -colorMapMax 1 -colorMapMin 0 -algorithm naview -sequenceDBN $seq -structureDBN $struct -colorMap $entropy -o ${newname}.new.svg


java -Djava.awt.headless=true -cp $now_dir/software_needed/bin/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd -algorithm naview -sequenceDBN $seq -structureDBN $struct -colorMap $entropy -o ${newname}.new.svg