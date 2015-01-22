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
if(length(commandArgs(TRUE))){
    chr <- commandArgs(TRUE)[1]
    verbose=FALSE
}
########################################################################
## Details
root <- "~/selection/counts/all"
out <- "~/selection/analysis/gscan/"
results.tag <- "read"
read.root <- "~/data/v6/reads/jj2"
indfile <- "~/data/v6/use/v61kg_europe2names.ind"

pops <- c("WHG", "EN", "Yamnaya", "CEU", "GBR", "IBS", "TSI")
A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.226, 0, 0.712, 0.287),3, 4) 

########################################################################

include.reads <- list(                  #Include these populations as reads
    "WHG"=c("SpanishMesolithic", "HungaryGamba_HG"), #SpanishMesolithic is the high coverage LaBrana I0585
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

## get list of samples in each population of reads
include.read.samples <- lapply(include.reads, function(x){NULL})
for(pop in names(include.reads)){
    for(subpop in include.reads[[pop]]){
        include.read.samples[[pop]] <- c(include.read.samples[[pop]], ind[ind[,3]==subpop,1])
    }
}

## Sort reads by ID for faster indexing
reads <- reads[order(match(reads[,1], data$ID)),]
read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
N.read.samples <- as.numeric(read.sample.counts[1])
if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}


## set up results
results <- matrix(0, nrow=NROW(data), ncol=2)
rownames(results) <- data$ID

empty.data <- rep( list(list()), length(pops) ) 
names(empty.data) <- pops
for(pop in pops){
    empty.data[[pop]] <- list("reads"=list("ref"=NULL, "alt"=NULL), "counts"=c(0,0)) #ref and alt counts. 
}


for(i in 1:NROW(data)){
    this.snp <- data[i,1]
    if(verbose){cat(paste0("\r", i, " ", this.snp))}

    ## this.read <- reads[reads[,1]==this.snp,]
    ## Select the read counts for this snp, but double check that we've got the right ones!
    this.read <- reads[(1+N.read.samples*(i-1)):(N.read.samples*i),]
    if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

                                        #Setup read data
    freq.data <- empty.data
    for(pop in pops){
        if(pop %in% names(include.reads)){
            for(sample in include.read.samples[[pop]]){
                ref.alt <- this.read[this.read[,2]==sample,3:4]
                if(sum(ref.alt)>0){
                    freq.data[[pop]][["reads"]][["ref"]] <- c(freq.data[[pop]][["reads"]][["ref"]],ref.alt[[1]])
                    freq.data[[pop]][["reads"]][["alt"]] <- c(freq.data[[pop]][["reads"]][["alt"]],ref.alt[[2]])
                }
            }
        }

        if(pop %in% names(include.counts)){
            for(subpop in include.counts[[pop]]){
                ref.alt <- c(counts[i,subpop], totals[i,subpop]-counts[i,subpop])
                names(ref.alt) <- NULL
                freq.data[[pop]][["counts"]] <- freq.data[[pop]][["counts"]]+ref.alt
            }
        }
    }    
    
    results[i,] <- test.3pop.reads(freq.data, A)
}

results <- cbind(rownames(results), results)
colnames(results) <- c("ID", "ChiSq", "uncorrected.p")
results <- data.frame(results)
write.table(results, paste0("~/selection/analysis/gscan/scan_results_read", results.tag, ".chr", chr, ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

