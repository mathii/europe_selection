#This is the genome-wide scan but using the read level information to try and
#get some idea about diploid calls. For reasons of speed, we break this one
#by chromosome. Genomic correction has to be done in another step for this reason.

#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
source("~/selection/code/lib/3pop_lib.R")

#Modern GBR, CEU, IBS, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################

source("~/selection/code/analysis/setup_populations_reads.R")

########################################################################

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
results <- matrix(0, nrow=NROW(data), ncol=5)
rownames(results) <- data$ID

## Data structure
empty.data <- make.empty.data(pops)

## Check monocheck
if(!all(monocheck %in% colnames(counts))){
    stop(paste0("Missing populations: ", monocheck[!(monocheck %in% colnames(counts))]))
}

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
        eff.N <-  round(effective.data.size(freq.data)[1:3], 2)
        results[i,] <- c(test.3pop.reads(freq.data, A, error.prob=error.prob), eff.N)
    }
}

results <- results[!is.na(results[,2]),]

results <- cbind(rownames(results), results)
colnames(results) <- c("ID", "ChiSq", "uncorrected.p", "eff.N1", "eff.N2", "eff.N3")
results <- data.frame(results)
out.file <-  paste0("~/selection/analysis/",version,"/gscan/scan_results_read", results.tag, ".chr", chr, ".txt")
print(out.file)
write.table(results,out.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

