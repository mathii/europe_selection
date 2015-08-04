## Make a table summarising the LCT imputation results.

results <- read.table("~/selection/analysis/v8/LCT/rs4988235.imputation_results", as.is=TRUE, header=TRUE, comment="")
results <- results[,10:NCOL(results)]
inc <- grepl("0|1", results[2,], fixed=TRUE)|grepl("1|0", results[2,], fixed=TRUE)|grepl("1|1", results[2,], fixed=TRUE)
results <- results[,inc]

imputed.samples <- names(results)
imputed.samples[which(imputed.samples=="BR2")] <- "I1504"

## Load reads and data
reads <- read.table(paste0("~/data/v8/reads/jj2.chr2.readcounts.gz"), as.is=TRUE, header=FALSE)
ind <- read.table("~/data/v8/use/v81kg_europe2names.ind", as.is=TRUE, header=FALSE)
rownames(ind) <- ind[,1]

tab <- data.frame(sample=imputed.samples, stringsAsFactors=FALSE)
reads <- reads[reads[,1]=="rs4988235",]
rownames(reads) <- reads[,2]
tab$pop <- ind[tab$sample,3]
tab$ref <- reads[tab$sample,3]
tab$alt <- reads[tab$sample,4]

