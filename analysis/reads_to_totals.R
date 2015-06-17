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

#Include these populations as hard calls
include.totals <- c( "Loschbour", "Stuttgart", "CEU", "GBR", "IBS", "TSI", "YRI")

## Setup the data. 
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
totals <- totals[totals[,"CHR"]%in%chrs,]
totals <- totals[order(totals$CHR),]
data <- totals[,1:5]
totals <- totals[,6:NCOL(totals)]

new.totals <- 0*totals
## new.totals$SpanishMesolithic <- 0

totals <- data.matrix(totals)
new.totals <- data.matrix(new.totals)

## get list of samples in each population of reads
ind <- read.table(indfile, as.is=TRUE, header=FALSE)
include.read.samples <- list()
for(j in 1:NCOL(new.totals)){
    if( !(colnames(new.totals)[j] %in% include.totals)){
        pop <- colnames(new.totals)[j]
        include.read.samples[[pop]] <- ind[ind[,3]==pop,1]
    }
}

this.chr=0                             #Which chromosome are we currently on?
readi <- 1
for(i in 1:NROW(data)){
    if(this.chr!=data[i,"CHR"]){
        this.chr <- data[i,"CHR"]
        cat(paste0("Loading chromosome ", this.chr, " reads..."))
        reads <- read.table(paste0(read.root, ".chr", data[i,"CHR"], ".readcounts.gz"), as.is=TRUE, header=FALSE)
        #Restrict to SNPs that are in data
        reads<-reads[reads[,1] %in% data[,1],]
        this.chr.ID.order <- data[data[,"CHR"]==this.chr,"ID"]
        reads <- reads[order(match(reads[,1],this.chr.ID.order)),]
        read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
        N.read.samples <- as.numeric(read.sample.counts[1])
        if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}
        readi <- 1
        cat(paste0("Done\n"))
    }

    
    if(verbose){cat(paste0("\r", readi))}
    this.snp <- data[readi,1]
    ## this.read <- reads[reads[,1]==this.snp,]
    this.read <- reads[(1+N.read.samples*(readi-1)):(N.read.samples*readi),]
    if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

    for(j in 1:NCOL(new.totals)){
        if(colnames(new.totals)[j] %in% include.totals){
            new.totals[i,j] <- totals[i,j]
        } else{
            pop=colnames(new.totals)[j]
            tt <- 0 
            for(sample in include.read.samples[[pop]]){
                if(any(this.read[,2]==sample)){
                    ref.alt <- this.read[this.read[,2]==sample,3:4]
                    tt <- tt+2-0.5^(ref.alt[1]+ref.alt[2]-1)
                }
            }
            new.totals[i,j] <- tt
        } 
    }
    readi <- readi+1
}

write.table(new.totals, paste0("~/selection/analysis/", version, "/effsize/effsize_reads", ".chr", paste(chrs, collapse="_"), ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
