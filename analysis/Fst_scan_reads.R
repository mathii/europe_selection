#Genome-wide scan for differences in Fst
# First calculate the allele frequencies based on the reads. 

source("~/selection/code/lib/3pop_lib.R")

########################################################################
## Details
chr <- 1                                #set manually, or from --args
verbose=TRUE
version <- "vx" #v6, v7 etc...

cA <- commandArgs(TRUE)
if(length(cA)){
  chr <- cA[1]
  version <- cA[2]
  ## verbose=FALSE
}

########################################################################
## Details
root <- paste0("~/selection/counts/", version, "/all")
out <- paste0("~/selection/analysis/", version, "/gscan/")
read.root <- paste0("~/data/", version, "/reads/jj2")
indfile <- paste0("~/data/", version, "/use/", version,"1kg_europe2names.ind")
error.prob <- 0.001

pops <- c("HG", "EN", "MN", "LNBA", "SBA", "CEU", "GBR", "IBS", "TSI", "FIN")
#Check if the SNP is monomorphic in these populations. 
monocheck <- c("CEU", "GBR", "IBS", "TSI", "HungaryGamba_HG", "Loschbour", "Stuttgart",
               "LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN", "Yamnaya")

########################################################################

include.reads <- list(                  #Include these populations as reads
    "HG"=c("SpanishMesolithic", "Motala_HG", "HungaryGamba_HG"), 
    "EN"=c("LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN"), 
    "MN"=c("Esperstedt_MN", "Spain_MN", "Iceman", "Baalberge_MN"),
    "LNBA"=c("Bell_Beaker_LN", "BenzigerodeHeimburg_LN", "Corded_Ware_LN", "Halberstadt_LBA", "HungaryGamba_BA", "HungaryGamba_CA"),
                      "SBA"=c("Yamnaya", "Poltavka", "Srubnaya")
                      )
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "EN"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI", "FIN"="FIN" )

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
reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
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
results <- matrix(0, nrow=NROW(data), ncol=3)
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
        f <- fit.unconstrained.model.reads(freq.data)$par
        ## w <- effective.data.size(freq.data)
        ## w <- w/sum(w)
        Hs <- mean(2*f*(1-f))
        Ht <- 2*mean(f)*(1-mean(f))
        results[i,] <- c(Hs, Ht, 1-Hs/Ht)
    }
}

results <- results[!is.na(results[,2]),]

results <- cbind(rownames(results), results)
colnames(results) <- c("ID", "Hs", "Ht", "Fst")
results <- data.frame(results)

out.file <-  paste0("~/selection/analysis/",version,"/fst/fst_scan_results_read.chr", chr, ".txt")
print(out.file)
write.table(results,out.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

