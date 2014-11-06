anno <- read.table("V4_transcript_annotation_0331.txt", sep = '\t', header=T, stringsAsFactors=F)
bioType <- sort(table(anno[['gene_biotype']]), decreasing=T)
src <- sort(table(anno[['source']]), decreasing=T)
write.table(bioType, file = 'bio_type.txt', col.names=F, quote=F, sep = '\t')
write.table(src, file = 'source.txt', col.names=F, quote=F, sep = '\t')

anno[anno[, 3] == '-', c(1, 3, 4, 5)]

anno[anno[, 3] == '-', 3] <- 'lincRNA'

write.table(anno[, c(1,2,3,6,7)], file = 'V4_final_transcript_annotation_0401.txt', 
            sep = '\t', quote = F, row.names = F)

final <- anno[, c(1,2,3,6,7)]
