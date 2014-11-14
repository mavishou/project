#!/bin/bash 

# -------------get the final list and sequence of lncrnadb lncRNAs that are not existed in GENCODE v19-------------------

# get lncRNAs that could be mapped to lncrnadb in v21 but not exist in v19 
myjoin all_trans.txt  mapped_gencode_lncrnas.txt  | grep "^-" | cut -f 2 > not_in_v19

# get lncRNAs that could be mapped to lncrnadb in v21 and exist in v19
myjoin all_trans.txt  mapped_gencode_lncrnas.txt  | grep "^=" | cut -f 2 > in_v19

# since transcripts not_in_v19 has version, it could the reason that not in in v19, so remove the version
cat not_in_v19 | sed -r 's/\.[0-9]+//' > not_in_v19_no_version

# get all_trans no version
cat all_trans.txt | sed -r 's/\.[0-9]+//' > all_trans_no_version

# get gencode v21 lncRNAs that mapped to lncrnadb but could not be found in v19, in spite of version
# only one lncRNA found: ENST00000613780
myjoin all_trans_no_version not_in_v19_no_version | grep "^-" | cut -f 2

# which lncrnadb lncRNA does it influence?
grep ENST00000613780 gencode_whole_lncRNAdb.final.format
# the reuslt:
#ENST00000613780.1       gomafu_homosapiens_1

# then, modify the not aligned lncrnadb lncrnas sequence file, and add the gomafu_homosapiens_1 to the unmapped list
cp js/not.aligned.lncRNAdb.fa .
mv not.aligned.lncRNAdb.fa not_aligned_lncRNAdb.fa
vim not_aligned_lncRNAdb.fa
# remove linc00237_hg_1 since we have no sequence
# add gomafu_homosapiens_1

#-------------------------------------------test blat--------------------------------------

# start the blat gfserver
gfServer start blatMachine 33333 -tileSize=11 -log=blatServerCi1.log /rd/user/houm/genome/human/hg19/hg19.2bit

# begin blat
gfClient -t=dna -q=rna 127.0.0.1 33333 / test.fa test.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl test.psl out.psr
perl psl2bed.pl test.psl test.all.bed12
bedtools bed12tobed6 -i test.all.bed12 > test.bed

# test for single exon
gfClient -t=dna -q=rna 127.0.0.1 33333 / se.fa se.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl se.psl out.psr
perl psl2bed.pl se.psl se.all.bed12
bedtools bed12tobed6 -i se.all.bed12 > se.bed
# result: the exon end is plused by 1, just as what TX said

# test for multi exon
gfClient -t=dna -q=rna 127.0.0.1 33333 / me.fa me.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl me.psl out.psr
perl psl2bed.pl me.psl me.all.bed12
bedtools bed12tobed6 -i me.all.bed12 > me.bed
# result: correct. Pay attention: if there are several best matches, the result will include them all.

# test for mismatch
gfClient -t=dna -q=rna 127.0.0.1 33333 / 1mm.fa 1mm.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl 1mm.psl out.psr
perl psl2bed.pl 1mm.psl 1mm.all.bed12
bedtools bed12tobed6 -i 1mm.all.bed12 > 1mm.bed
# result: The bed is the same as no mismatch

# test for insertion in query seq
gfClient -t=dna -q=rna 127.0.0.1 33333 / is.fa is.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl is.psl out.psr
perl psl2bed.pl is.psl is.all.bed12
bedtools bed12tobed6 -i is.all.bed12 > is.bed
# result: the target seq (genome) is splitted into 2 parts:
# chr1    69090   69270   ENST00000335137 917     +
# chr1    69270   70008   ENST00000335137 917     +

# test for small case sequence
gfClient -t=dna -q=rna 127.0.0.1 33333 / sc.fa sc.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl sc.psl out.psr
perl psl2bed.pl sc.psl sc.all.bed12
bedtools bed12tobed6 -i sc.all.bed12 > sc.bed

# ===================get the gtf of lnrnadb lncRNAs not in GECODE v19========================
#------------------------run blat for lncrnadb sequence -------------------
cd /rd/user/houm/projects/AnnoLnc/gene_model/blat
gfClient -t=dna -q=rna 127.0.0.1 33333 / ../not_aligned_lncRNAdb.fa not_aligned_lncRNAdb.psl
pslSort dirs sorted.psl temp ./
pslReps sorted.psl final.psl out.psr
perl /rd/user/houm/projects/AnnoLnc/bin/psl2bed.pl final.psl final.all.bed12

# ------------analyze and check the blat result--------------------
# how many not aligned lncRNAs in total? 74
 grep "^>" ../not_aligned_lncRNAdb.fa | wc -l

# how many lncRNAs are mapped to the human genome? 74
cut -f 4 final.all.bed12 | sort | uniq | wc -l









