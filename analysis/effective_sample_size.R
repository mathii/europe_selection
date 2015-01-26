#Calculate the effective sample size
source("~/selection/code/lib/readlib.R")

########################################################################
## Details
root <- "~/selection/counts/all"
read.totals <- "~/selection/analysis/effsize/effsize_total.txt"
########################################################################

## Majority call info
rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

## Read level info
read.data <- read.table(read.totals, as.is=TRUE, header=TRUE)
read.totals <- read.data[,6:NCOL(read.data)]
read.data <- read.data[,1:5]



results <- data.frame(pops=colnames(totals), size=apply(totals, 2, max), effective.size=colMeans(totals))
write.table(results, "~/selection/analysis/effective_sample_size.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

