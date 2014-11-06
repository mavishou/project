setwd('/lustre/user/houm/projects/AnnoLnc/CLIP/pipeclip/analysis/0.01')
intersectFile <- read.table('intersect_FUS.txt', sep='\t', stringsAsFactors=F)

# make the cluster annotation
clusterAnno <- intersectFile[, 1:9]
clusterAnno <- unique(clusterAnno)
clusterLen <- clusterAnno[, 3] - clusterAnno[, 2]
clusterReadsNum <- clusterAnno[, 5]
clusterFDR <- clusterAnno[, 7]
names(clusterLen) <- names(clusterReadsNum) <- names(clusterFDR) <- clusterAnno[, 4]

cluster2trans <- cbind(intersectFile[, 4], 
                       do.call('rbind', strsplit(intersectFile[, 18], ' '))[, 4])

# make the transcript annotation
transAnno <- read.table('/lustre/user/houm/projects/AnnoLnc/V4_final_transcript_annotation_0401.txt', sep='\t', stringsAsFactors=F, header=T)
if(any(! cluster2trans[, 2] %in% transAnno[, 1])){
  stop('transcript not existed!')
}
transGeneId <- transAnno[, 2]
transType <- transAnno[, 3]
transGeneName <- transAnno[, 4]
names(transGeneId) <- names(transType) <- names(transGeneName) <- transAnno[, 1]

rm(intersectFile, clusterAnno, transAnno)

outClusterAnno <- 
  cbind(tapply(cluster2trans[, 1], as.factor(cluster2trans[, 2]), length), 
        tapply(cluster2trans[, 1], as.factor(cluster2trans[, 2]), function(x){
          paste(x, collapse='; ') # cluster name
        }),
        tapply(cluster2trans[, 1], as.factor(cluster2trans[, 2]), function(x){
          paste(clusterReadsNum[x], collapse='; ') # cluster reads num
        }),
        tapply(cluster2trans[, 1], as.factor(cluster2trans[, 2]), function(x){
          sum(clusterReadsNum[x]) # cluster reads num sum
        }),
        tapply(cluster2trans[, 1], as.factor(cluster2trans[, 2]), function(x){
          paste(clusterFDR[x], collapse='; ') # cluster FDR
        }),
        tapply(cluster2trans[, 1], as.factor(cluster2trans[, 2]), function(x){
          min(clusterFDR[x]) # min cluster FDR
        })
  )


outTrans <- rownames(outClusterAnno)
outAnno <- cbind(transGeneName[outTrans], outTrans, transGeneId[outTrans], outClusterAnno[, c(1, 4, 6, 2, 3, 5)])
# outAnno <- cbind(transGeneName[outTrans], outTrans, transGeneId[], outClusterAnno)
colnames(outAnno) <- c('GeneName', 'TranscriptID', 'GeneID', 'ClusterNum', 
                       'SumClusterReadsNum', 'MinClusterFDR', 'ClusterNames', 
                       'ClusterReadsNum', 'ClusterFDRs')

outAnno <- outAnno[order(-as.numeric(outAnno[, 'SumClusterReadsNum']), as.numeric(outAnno[, 'MinClusterFDR'])), ]
outPath <- ''
write.table(outAnno, file=outPath, sep='\t', quote=F, row.names=F)

# stat
numCluster <- length(clusterLen)
numTrans <- length(outTrans)
numGenes <- length(unique(outAnno[, 'GeneID']))

cat('# clusters that are annotated: ')
cat (numCluster)
cat('\n')
cat('# Transcripts: ')
cat (numTrans)
cat('\n')
cat('# Genes: ')
cat (numGenes)
cat('\n')

