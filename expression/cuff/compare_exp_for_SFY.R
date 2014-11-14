# read the FPKM table
fpkmTableHM <- read.table("/lustre/user/houm/projects/AnnoLnc/expression/cuff/cuffnorm/isoforms.fpkm_table", sep = "\t", header = T, stringsAsFactors=F)
fpkmTableSFY <- read.table("/lustre/user/shify/expression/bam/norm_out/isoforms.fpkm_table", sep = "\t", header = T, stringsAsFactors=F)

# check
all(fpkmTableHM[, 1] == fpkmTableSFY[, 1])

# select sample
hm <- fpkmTableHM[, 2]
sfy <- fpkmTableSFY[, 3]

# calculate pearson correlation coeffecient
cor(hm, sfy)

# plot
plot(hm, sfy, pch=20)
plot(hm, sfy, xlim=c(0, 1000), ylim=c(0, 1000), pch=20, 
     xlab="HM", ylab="SFY", 
     main="Compare ERR030882 transcript FPKM calculated by SFY and HM")
