## Merge imputed and non-imputed results

cA <- commandArgs(TRUE)
imputed <- cA[1]
reads <- cA[2]
out <- cA[3]
cutoff <- as.numeric(cA[4])

imputed <- read.table(imputed, header=TRUE, as.is=TRUE)
reads <- read.table(reads, header=TRUE, as.is=TRUE)

mrg<-merge(imputed, reads, by="ID", ALL=TRUE, suffixes=c(".i", ".r"))
mrg$ChiSq <- ifelse(mrg$AR2>cutoff, mrg$ChiSq.i, mrg$ChiSq.r)
mrg$uncorrected.p <- ifelse(mrg$AR2>cutoff, mrg$uncorrected.p.i, mrg$uncorrected.p.r)
write.table(mrg, out, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
cat(paste0( "Included ", sum(mrg$AR2>cutoff), " imputed SNPs\n"))
