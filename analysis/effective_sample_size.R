#Calculate the effective sample size
source("~/selection/code/lib/readlib.R")

version=NA
if(length(commandArgs(TRUE))){
  version=commandArgs(TRUE)[1]
}

########################################################################
## Details
root <- paste0("~/selection/counts/", version ,"/all" )
read.totals <- paste0("~/selection/analysis/",version,"/effsize/effsize_reads.txt.gz")
########################################################################

## Majority call info
rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

## Read level info
read.totals <- read.table(read.totals, as.is=TRUE, header=TRUE)



results <- data.frame(pops=colnames(totals), size=apply(totals, 2, max), effective.size=colMeans(totals))
write.table(results, paste0("~/selection/analysis/",version,"/effsize/effective_sample_size.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

read.results <- data.frame(pops=colnames(read.totals), effective.size=colMeans(read.totals))
write.table(read.results, paste0("~/selection/analysis/",version,"/effsize/effective_sample_size_reads.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
