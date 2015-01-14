#Calculate the effective sample size
source("~/selection/code/lib/readlib.R")

########################################################################
## Details
root <- "~/selection/counts/all"

########################################################################

rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

results <- data.frame(pops=colnames(totals), size=apply(totals, 2, max), effective.size=colMeans(totals))
write.table(results, "~/selection/analysis/effective_sample_size.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
