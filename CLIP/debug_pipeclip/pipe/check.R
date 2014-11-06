cl <- read.table('2.filter.cluster.bed', sep='\t', stringsAsFactors=F)
l <- cl[, 3] - cl[, 2]
r <- cl[, 5]
ratio <- r/l
plot(ratio*100, cl[, 7]*100, xlim=c(0, 100))
