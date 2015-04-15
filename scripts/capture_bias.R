## Test snps for capture bias.

source("~/selection/code/lib/3pop_lib.R")

#Modern GBR, CEU, IBS, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
chr <- 22                                #set manually, or from --args
verbose=TRUE
if(length(commandArgs(TRUE))){
    chr <- commandArgs(TRUE)[1]
    verbose=FALSE
}

########################################################################
## Details
read.root <- "~/data/v6/reads/jj2"
error.prob <- 0.001
alpha <- 0.05
########################################################################

ll.diff <- qchisq(0.05, df=1, lower.tail=F)/2

reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts"), as.is=TRUE, header=FALSE)
## Must sort!
reads <- reads[order(reads[,1]),]

## Sort reads by ID for faster indexing
read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
N.read.samples <- as.numeric(read.sample.counts[1])
N.SNPS <- length(read.sample.counts)
if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}

results <- data.frame(SNP=rep("",N.SNPS), max.f=rep(0,N.SNPS), max.p=rep(0,N.SNPS), low.p=rep(0,N.SNPS), hi.p=rep(0, N.SNPS), stringsAsFactors=FALSE)
empty.data <- make.empty.data("ALL")

for(i in 1:N.SNPS){
    if(verbose){cat(paste0("\r", i, " ", this.read[1,1]))}

    this.read <- reads[(1+N.read.samples*(i-1)):(N.read.samples*i),]
    if(!all(this.read[,1]==this.read[1,1])){stop("Selected the wrong SNP")}
    
    data <- empty.data
    data[["ALL"]][["reads"]][["ref"]]=this.read[,3]
    data[["ALL"]][["reads"]][["alt"]]=this.read[,4]

    these.log.likelihoods <- matrix(0, nrow=11, ncol=101)
    for(k in 1:11){
        for(j in 1:101){
            hp=min(max((j-1)/100, error.prob), 1-error.prob) 
            these.log.likelihoods[k,j] <- likelihood.reads((k-1)/10, data, error.prob=error.prob, het.p=hp)
        }
    }
    
    best.ij <- which(these.log.likelihoods==max(these.log.likelihoods), arr.ind=TRUE)[1,]
    max.ll <- these.log.likelihoods[best.ij[1], best.ij[2]]
    low.pi=best.ij[2]
    hi.pi=best.ij[2]
    for(k in 1:11){
        if(any(these.log.likelihoods[k,]>max.ll-ll.diff)){
            low.pi=min(low.pi, min(which(these.log.likelihoods[k,]>max.ll-ll.diff)))
            hi.pi=max(hi.pi, max(which(these.log.likelihoods[k,]>max.ll-ll.diff)))
        }
    }

    results[i,] <- c(this.read[1,1], (best.ij[1]-1)/10, (best.ij[2]-1)/100, (low.pi-1)/100, (hi.pi-1)/100)
}

write.table(results, paste0("~/capture_bias/390k_chr", chr, ".results"), row.names=FALSE, col.names=TRUE, quote=FALSE)

