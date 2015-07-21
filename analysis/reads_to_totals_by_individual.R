## Use the read count information to estimate the effective size of the data.
## Should really be implemented in spindrift. Actually this just counts
## 2 if we saw both alleles, and 1=0.5^n.reads if we only saw one. 
source("~/selection/code/lib/3pop_lib.R")

#Modern GBR, CEU, IBS, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
chrs <- 1:22                                #set manually, or from --args
verbose=TRUE
if(length(commandArgs(TRUE))){
    cc <- commandArgs(TRUE)[1]
    if(cc=="All"){
        chrs <- 1:22
    } else{
        chrs <- as.numeric(cc)
    }
    version <- commandArgs(TRUE)[2]
    verbose=FALSE
}


########################################################################
## Details
root <- paste0("~/selection/counts/", version, "/all")
out <- paste0("~/selection/counts/", version, "/all.reads")
read.root <- paste0("~/data/", version, "/reads/jj2")
indfile <- paste0("~/data/", version, "/use/", version, "1kg_europe2names.ind")

########################################################################

include.totals <- c( "Loschbour", "Stuttgart", "CEU", "GBR", "IBS", "TSI", "YRI", "FIN")

## Setup the data. 
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
totals <- totals[totals[,"CHR"]%in%chrs,]
totals <- totals[order(totals$CHR),]
data <- totals[,1:5]
totals <- totals[,6:NCOL(totals)]

new.totals <- 0*totals
## new.totals$SpanishMesolithic <- 0

totals <- data.matrix(totals)

## get list of samples in each population of reads
ind <- read.table(indfile, as.is=TRUE, header=FALSE)
include.read.samples <- c()
for(j in 1:NCOL(totals)){
    if( !(colnames(totals)[j] %in% include.totals)){
        pop=colnames(totals)[j]
        samples <-  ind[ind[,3]==pop,1]
        for(samp in samples){
            include.read.samples <- c(include.read.samples, samp)
        }
    }
}
include.read.samples <- sort(include.read.samples)

new.totals <- matrix(0, ncol=length(include.read.samples), nrow=nrow(data))
new.totals <- data.matrix(new.totals)
colnames(new.totals) <- include.read.samples

this.chr=0                             #Which chromosome are we currently on?
readi <- 1
for(i in 1:NROW(data)){
    if(this.chr!=data[i,"CHR"]){
        this.chr <- data[i,"CHR"]
        cat(paste0("Loading chromosome ", this.chr, " reads..."))
        reads <- read.table(paste0(read.root, ".chr", data[i,"CHR"], ".readcounts.gz"), as.is=TRUE, header=FALSE)
        this.chr.ID.order <- data[data[,"CHR"]==this.chr,"ID"]
        reads <- reads[order(match(reads[,1],this.chr.ID.order), reads[,2]),]
        reads <- reads[reads[,2]%in%include.read.samples,]
        read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
        N.read.samples <- as.numeric(read.sample.counts[1])
        if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}
        readi <- 1
        cat(paste0("Done\n"))
    }
    
    if(verbose){cat(paste0("\r", readi))}
    this.snp <- data[i,1]
    ## this.read <- reads[reads[,1]==this.snp,]
    this.read <- reads[(1+N.read.samples*(readi-1)):(N.read.samples*readi),]
    if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

    if(verbose){cat(paste0("\r", i))}
    new.totals[i,] <- this.read[,3]+this.read[,4]
    readi <- readi+1
}

results <- data.frame(IID=include.read.samples, coverage=colMeans(new.totals), hit.once=colMeans(new.totals>0), effective=colMeans(2-0.5^(new.totals-1)))

tag=paste0("~/selection/analysis/", version, "/effsize/effsize_reads_by_ind", ".chr", paste(chrs, collapse="_"), ".txt")
if(all(chrs==1:22)){
    tag <- paste0("~/selection/analysis/", version, "/effsize/effsize_reads_by_ind.txt")
}
write.table(results, tag, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## pdf("~/Desktop/Coverage_vs_Hits.pdf")
## plot(results$coverage, results$hit.once, ylab="Proportion of sites hit at least once", xlab="Mean coverage", pch=16, col="#377EBA", log="x")
## xpts <- 10^seq(-2, 2, length.out=100)
## ypts=1-exp(-xpts)
## lines(xpts, ypts, col="red")
## abline(h=0.95, col="grey", lty=2)
## abline(v=3, col="grey", lty=2)
## abline(v=8, col="grey", lty=2)
## dev.off()
