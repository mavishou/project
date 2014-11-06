# ============================process====================================
setwd("~/BaiduYun/AnnoLnc/CLIP/")
oldSample <- read.xlsx2("clip_GEO_samples.xlsx", sheetName = "all", stringsAsFactors = F)

rbpTable <- read.xlsx2("clip_GEO_samples.xlsx", sheetName = "human_check_RBPs", colIndex = 1:2, stringsAsFactors = F, header = F)
rbpMap <- rbpTable[, 2]
names(rbpMap) <- rbpTable[, 1]

all(oldSample[, 6] %in% rbpTable[, 1])
# [1] TRUE

newSample <- cbind(oldSample, rbpMap[oldSample[, 6]])

write.xlsx2(newSample, "clip_GEO_samples.xlsx", sheetName = "integrate", append = T, row.names=F)

# analysis select
select <- read.xlsx2("clip_GEO_samples.xlsx", sheetName = "select", stringsAsFactors = F)
rbp <- select[, 6]
length(unique(rbp))
# [1] 64
series <- select[, 2]
length(unique(series))
# [1] 37
