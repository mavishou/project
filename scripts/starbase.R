clipSamples <- read.table('clip_samples.txt', sep = '\t', stringsAsFactors = F)
# 有85个sample
allPros <- unique(clipSamples[, 3])
