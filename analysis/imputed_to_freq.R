## Use the read count information to estimate the freqeuncies in each population
## Should really be implemented in spindrift. Actually this just counts
## Here we're using the polymap populations. 
source("~/selection/code/lib/3pop_lib.R")
source("~/selection/code/lib/readlib.R")

########################################################################
## Details
chrs <- NA                                #set manually, or from --args
which.map <- ""
cA <- commandArgs(TRUE)
if(length(cA)){
    chr <- as.numeric(cA[1])
    version <- cA[2]
    if(length(cA)>2){
        which.map <- cA[3]
    }
}

verbose=TRUE
## Supposed to check if running on cluster, but YMMV
if( Sys.info()["login"]!=Sys.info()["user"]){
    verbose=FALSE
}

########################################################################
## Details
root <- paste0("~/selection/counts/",version,"/all")
impute.file <- paste0("~/selection/imputation/", version, "/imputed.within.chr", chr, ".vcf.gz")
indfile <- paste0("~/data/",version,"/use/v61kg_europe2names.ind")
snpfile <- paste0("~/data/",version,"/use/v61kg_europe2names.snp")
polymap <- "~/selection/code/files/polymap.txt"
out <- paste0("~/selection/counts/",version,"/all.reads")
if(which.map!=""){
    polymap <- paste0("~/selection/code/files/polymap.", which.map, ".txt" )
    out <- paste0("~/selection/counts/",version,"/all.reads.", which.map)
}

error.prob <- 0.01
########################################################################
#Include these populations as hard calls
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "Germany_EN"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
#Exclude these
exclude <- c("LaBrana1")
## exclude <- c()
                                        # include these
include.extra <- list("SpanishMesolithic"="WHG")         #High coverage LaBrana
monocheck <- c("CEU", "IBS", "GBR", "IBS", "TSI")
########################################################################

include.probs <- list()
polymap <- read.table(polymap, as.is=TRUE, header=FALSE)
for(i in 1:NROW(polymap)){
    if(polymap[i,2] %in% exclude){next}
    
    if(polymap[i,2] %in% names(include.probs)){
        include.probs[[polymap[i,2]]] <- c(include.probs[[polymap[i,2]]], polymap[i,1])
    } else{
        include.probs[[polymap[i,2]]] <- polymap[i,1]
    }
}

for(i in 1:length(include.extra)){
     if(include.extra[[i]] %in% names(include.probs)){
        include.probs[[include.extra[[i]]]] <- c(include.probs[[include.extra[[i]]]], names(include.extra)[i])
    } else{
        include.probs[[include.extra[[i]]]] <- names(include.extra)[i]
    }

}

pops <- unique(sort(c(names(include.probs), names(include.counts))))

## Setup the data. 
rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

include <- data$CHR == chr
counts <- counts[include,]
totals <- totals[include,]
data <- data[include,]

freq <- ci.low <- ci.up <- matrix(0, nrow=NROW(data), ncol=length(pops))
colnames(freq) <- colnames(ci.low) <- colnames(ci.up) <- pops
rownames(freq) <- rownames(ci.low) <- rownames(ci.up) <- data$ID

## get list of samples in each population of reads
include.prob.samples <- read.samples(indfile, include.probs, c(exclude, unlist(include.counts)))

########################################################################

#Load imputed likelihoods
impute <- read.table(impute.file, comment.char="", as.is=TRUE, header=FALSE, sep="\t", fill=TRUE)
comment.lines <- sum(grepl("^##", impute[,1]))
impute <- read.table(impute.file, comment.char="", as.is=TRUE, header=TRUE, skip=comment.lines, sep="\t", fill=TRUE)
rownames(impute) <- impute$ID
impute <- impute[impute$ID %in% data[,1],]
impute <- impute[,10:NCOL(impute)]

########################################################################


empty.data <- make.empty.data(pops)

for(i in 1:NROW(data)){
    if(verbose){cat(paste0("\r", i))}
    this.snp <- data[i,1]
    if(verbose){cat(paste0("\r", i, " ", this.snp))}
    this.prob <- impute[i,]

    freq.data <- make.prob.freq.data(pops, include.probs, include.prob.samples, include.counts,
                                this.prob, counts[i,], totals[i,], empty.data)

    monomorphic <- all(counts[i,monocheck]==0)|all(counts[i,monocheck]==totals[i,monocheck])
    if(monomorphic){
        freq[i,] <- NA
    }else{
        info <- ci.unconstrained.model.reads(freq.data, error.prob=error.prob)
        freq[i,] <- info$p
        ci.low[i,] <- info$lci
        ci.up[i,] <- info$uci
    }

}

freq <- cbind(data[,1:5],freq)
ci.low <- cbind(data[1:5],ci.low)
ci.up <- cbind(data[1:5],ci.up)
include <- !is.na(freq[,6])
freq <- freq[include,]
ci.low <- ci.low[include,]
ci.up <- ci.up[include,]

write.table(freq, paste0(out, ".chr", chr, ".imputed.freq"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(ci.low, paste0(out, ".chr", chr, ".imputed.lowCI.freq"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(ci.up, paste0(out, ".chr", chr, ".imputed.highCI.freq"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
