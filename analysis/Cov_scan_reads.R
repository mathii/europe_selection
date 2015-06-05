#Genome-wide scan for Covariance outliers. 
# First calculate the allele frequencies based on the reads. 

source("~/selection/code/lib/3pop_lib.R")

########################################################################
## Details
verbose=TRUE
version <- "vx" #v6, v7 etc...

cA <- commandArgs(TRUE)
if(length(cA)){
  version <- cA[1]
  ## verbose=FALSE
}

########################################################################
## Details
root <- paste0("~/selection/counts/", version, "/all")
out <- paste0("~/selection/analysis/", version, "/gscan/")
read.root <- paste0("~/data/", version, "/reads/jj2")
indfile <- paste0("~/data/", version, "/use/", version,"1kg_europe2names.ind")
error.prob <- 0.001

## TODO: Add FIN
pops <- c("HG", "EN", "MN", "LNBA", "SBA", "CEU", "GBR", "IBS", "TSI" )

########################################################################

degf <- length(pops)-1

## Setup the data. 
freqs <- read.table(paste0(root, ".reads.LargeScale.freq"), header=TRUE, as.is=TRUE)
freqs.hci <- read.table(paste0(root, ".reads.LargeScale.highCI.freq"), header=TRUE, as.is=TRUE)
freqs.lci <- read.table(paste0(root, ".reads.LargeScale.lowCI.freq"), header=TRUE, as.is=TRUE)

data <- freqs[,1:5]
rownames(freqs) <- data$ID
rownames(freqs.lci) <- freqs.lci$ID
freqs.lci <- freqs.lci[rownames(freqs),]
rownames(freqs.hci) <- freqs.hci$ID
freqs.hci <- freqs.hci[rownames(freqs),]


freqs <- data.matrix(freqs[,6:NCOL(freqs)])
freqs <- freqs[,pops]

freqs.hci <- data.matrix(freqs.hci[,6:NCOL(freqs.hci)])
freqs.hci <- freqs.hci[,pops]

freqs.lci <- data.matrix(freqs.lci[,6:NCOL(freqs.lci)])
freqs.lci <- freqs.lci[,pops]

include <- rowMeans(freqs)>0.05 & rowMeans(freqs)<0.95 & rowSums(1-freqs>1e-3)>2 & rowSums(freqs>1e-3)>2
include2 <- rowSums(freqs.hci-freqs.lci>0.55)==0
include <- include & include2
data <- data[include,]
freqs <- freqs[include,]
freqs.hci <- freqs.hci[include,]
freqs.lci <- freqs.lci[include,]

## Compute covariance on genome-wide SNPS?
ps <- freqs[,pops]
mn <- rowMeans(ps)
ps <- (ps-mn)/mn/(1-mn)
cov<-cov(ps)

K <- length(pops)
ps.drop <- ps[,1:(K-1)]
cov.inv.drop <- solve(cov[1:(K-1),1:(K-1)])
results <- matrix(0, nrow=NROW(data), ncol=2)
rownames(results) <- data$ID

for(i in 1:NROW(ps.drop)){
    cat(paste0("\r", i))
    stat <- ps.drop[i,,drop=F] %*% cov.inv.drop %*% t(ps.drop[i,,drop=F])
    results[i,] <- c(stat, pchisq(stat, df=degf, lower.tail=FALSE))
}
    

results <- results[!is.na(results[,2]),]

save.results <- cbind(rownames(results), results)

colnames(save.results) <- c("ID", "ChiSq", "uncorrected.p")
save.results <- data.frame(save.results)
out.file <-  paste0("~/selection/analysis/",version,"/cov/scan_results_reads.txt")
print(out.file)
write.table(save.results,out.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


## reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
## #Restrict reads to snps included in data file. 
## reads<-reads[reads[,1] %in% data[,1],]

## ## get list of samples in each population of reads
## include.read.samples <- read.samples(indfile, include.reads)

## ## Sort reads by ID for faster indexing
## reads <- reads[order(match(reads[,1], data$ID)),]
## read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
## N.read.samples <- as.numeric(read.sample.counts[1])
## if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}

## ## set up results
## results <- matrix(0, nrow=NROW(data), ncol=2)
## rownames(results) <- data$ID

## ## Data structure
## empty.data <- make.empty.data(pops)



## for(i in 1:NROW(data)){
##     this.snp <- data[i,1]
##     if(verbose){cat(paste0("\r", i, " ", this.snp))}

##     ## this.read <- reads[reads[,1]==this.snp,]
##     ## Select the read counts for this snp, but double check that we've got the right ones!
##     this.read <- reads[(1+N.read.samples*(i-1)):(N.read.samples*i),]
##     if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

##     #Setup read data
##     freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts,
##                                 this.read, counts[i,], totals[i,], empty.data)
##     monomorphic <- all(counts[i,monocheck]==0)|all(counts[i,monocheck]==totals[i,monocheck])
##     if(monomorphic){
##         results[i,] <- NA
##     }else{
##         f <- fit.unconstrained.model.reads(freq.data)$par
##         x=matrix(f-mean(f))
##         stat <- t(x) %*% solve(cov) %*% x
##         results[i,] <- c(stat, pchisq(stat, df=degf, lower.tail=FALSE))
        
##     }
## }

