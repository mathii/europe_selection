#This is the genome-wide scan but using the read level information to try and
#get some idea about diploid calls. For reasons of speed, we break this one
#by chromosome. Genomic correction has to be done in another step for this reason.

#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
source("~/selection/code/lib/3pop_lib.R")

#Modern GBR, CEU, IBS, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
chr <- 1                                #set manually, or from --args
verbose=TRUE
version <- "vx" #v6, v7 etc...
if(length(commandArgs(TRUE))){
    chr <- commandArgs(TRUE)[1]
    version <- commandArgs(TRUE)[2]
    verbose=FALSE
}
########################################################################
## Details
root <- paste0("~/selection/counts/", version, "all")
out <- paste0("~/selection/analysis/", version, "gscan/")
results.tag <- ""
read.root <- paste0("~/data/", version, "/reads/jj2")
indfile <- paste0("~/data/", version, "/use/", version,"1kg_europe2names.ind")
error.prob <- 0.001

pops <- c("WHG", "EN", "Yamnaya", "CEU", "GBR", "IBS", "TSI")
#Check if the SNP is monomorphic in these populations. 
monocheck <- c("CEU", "GBR", "IBS", "TSI", "LaBrana1", "HungaryGamba_HG", "Loschbour", "Stuttgart",
               "LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN", "Yamnaya")
A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.227, 0, 0.712, 0.288),3, 4) 

########################################################################

include.reads <- list(                  #Include these populations as reads
    "WHG"=c("LaBrana1", "HungaryGamba_HG"), #Replace LaBrana1 with SpanishMesolithic for the high coverage LaBrana I0585
    "EN"=c("LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN"), 
    "Yamnaya"="Yamnaya")
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "EN"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
## Setup the data. 
counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
include <- data$CHR==chr
data <- data[include,]

counts <- data.matrix(counts[,6:NCOL(counts)])
totals <- data.matrix(totals[,6:NCOL(totals)])
counts <- counts[include,]
totals <- totals[include,]
reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts"), as.is=TRUE, header=FALSE)
#Restrict reads to snps included in data file. 
reads<-reads[reads[,1] %in% data[,1],]

## get list of samples in each population of reads
include.read.samples <- read.samples(indfile, include.reads)

## Sort reads by ID for faster indexing
reads <- reads[order(match(reads[,1], data$ID)),]
read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
N.read.samples <- as.numeric(read.sample.counts[1])
if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}


## set up results
results <- matrix(0, nrow=NROW(data), ncol=2)
rownames(results) <- data$ID

## Data structure
empty.data <- make.empty.data(pops)

for(i in 1:NROW(data)){
    this.snp <- data[i,1]
    if(verbose){cat(paste0("\r", i, " ", this.snp))}

    ## this.read <- reads[reads[,1]==this.snp,]
    ## Select the read counts for this snp, but double check that we've got the right ones!
    this.read <- reads[(1+N.read.samples*(i-1)):(N.read.samples*i),]
    if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

                                        #Setup read data

    freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts,
                                this.read, counts[i,], totals[i,], empty.data)
    monomorphic <- all(counts[i,monocheck]==0)|all(counts[i,monocheck]==totals[i,monocheck])
    if(monomorphic){
        results[i,] <- NA
    }else{
        results[i,] <- test.3pop.reads(freq.data, A, error.prob=error.prob)
    }
}

results <- results[!is.na(results[,2]),]

results <- cbind(rownames(results), results)
colnames(results) <- c("ID", "ChiSq", "uncorrected.p")
results <- data.frame(results)
write.table(results, paste0("~/selection/analysis/",version,"/gscan/scan_results_read", results.tag, ".chr", chr, ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

