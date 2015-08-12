## Use the read count information to estimate the freqeuncies in each population
## Should really be implemented in spindrift. Actually this just counts
## Here we're using the polymap populations. 
source("~/selection/code/lib/3pop_lib.R")
source("~/selection/code/lib/readlib.R")

########################################################################
## Details
chrs <- 1:22                                #set manually, or from --args
which.map <- ""
cA <- commandArgs(TRUE)
if(length(cA)){
    chrs <- as.numeric(cA[1])
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
read.root <- paste0("~/data/",version,"/reads/jj2")
indfile <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.ind")
snpfile <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.snp")
polymap <- paste0("~/selection/code/files/",version, "/polymap.txt")
out <- paste0("~/selection/counts/",version,"/all.reads")
if(which.map!=""){
    polymap <- paste0("~/selection/code/files/", version, "/polymap.", which.map, ".txt" )
    out <- paste0("~/selection/counts/",version,"/all.reads.", which.map)
}

error.prob <- 0.01
########################################################################
#Include these populations as hard calls
always.include.counts <- c("Loschbour", "Stuttgart")
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "Germany_EN"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
#Exclude these
exclude <- c("LaBrana1")
## exclude <- c()
                                        # include these
if(version=="v6"){
include.extra <- list("SpanishMesolithic"="WHG")         #High coverage LaBrana
}else{
    include.extra <- list()
}

if(version=="v8"){
include.counts <- list(                 #Include these populations as hard calls. 
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
}

monocheck <- c("CEU", "IBS", "GBR", "FIN", "TSI", "YRI")
########################################################################

include.reads <- list()
polymap <- read.table(polymap, as.is=TRUE, header=FALSE)
for(i in 1:NROW(polymap)){
    if(polymap[i,2] %in% exclude){next}

    if( polymap[i,1] %in% always.include.counts ){
        if(polymap[i,2] %in% names(include.counts)){
            include.counts[[polymap[i,2]]] <- c(include.counts[[polymap[i,2]]], polymap[i,1])
        } else{
            include.counts[[polymap[i,2]]] <- polymap[i,1]
        }
    } else if(polymap[i,2] %in% names(include.reads)){
        include.reads[[polymap[i,2]]] <- c(include.reads[[polymap[i,2]]], polymap[i,1])
    } else{
        include.reads[[polymap[i,2]]] <- polymap[i,1]
    }
}

if(length(include.extra)){
    for(i in 1:length(include.extra)){
        if(include.extra[[i]] %in% names(include.reads)){
            include.reads[[include.extra[[i]]]] <- c(include.reads[[include.extra[[i]]]], names(include.extra)[i])
        } else{
            include.reads[[include.extra[[i]]]] <- names(include.extra)[i]
        }
    }
}

pops <- unique(sort(c(names(include.reads), names(include.counts))))

## Setup the data. 
rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

include <- data$CHR %in% chrs
counts <- counts[include,]
totals <- totals[include,]
data <- data[include,]

freq <- ci.low <- ci.up <- matrix(0, nrow=NROW(data), ncol=length(pops))
colnames(freq) <- colnames(ci.low) <- colnames(ci.up) <- pops
rownames(freq) <- rownames(ci.low) <- rownames(ci.up) <- data$ID

## get list of samples in each population of reads
include.read.samples <- read.samples(indfile, include.reads, c(exclude, unlist(include.counts)))

########################################################################

empty.data <- make.empty.data(pops)

this.chr=0                             #Which chromosome are we currently on?
readi <- 1
for(i in 1:NROW(data)){
    if(this.chr!=data[i,"CHR"]){
        this.chr <- data[i,"CHR"]
        cat(paste0("Loading chromosome ", this.chr, " reads..."))
        reads <- read.table(paste0(read.root, ".chr", data[i,"CHR"], ".readcounts.gz"), as.is=TRUE, header=FALSE)
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

    freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts,
                                this.read, counts[i,], totals[i,], empty.data)

    monomorphic <- all(counts[i,]==0)|all(counts[i,monocheck]==totals[i,monocheck])
    if(monomorphic){
        freq[i,] <- NA
    }else{
        info <- ci.unconstrained.model.reads(freq.data, error.prob=error.prob)
        freq[i,] <- info$p
        ci.low[i,] <- info$lci
        ci.up[i,] <- info$uci
    }

    readi <- readi+1
}

freq <- cbind(data[,1:5],freq)
ci.low <- cbind(data[1:5],ci.low)
ci.up <- cbind(data[1:5],ci.up)
include <- !is.na(freq[,6])
freq <- freq[include,]
ci.low <- ci.low[include,]
ci.up <- ci.up[include,]

write.table(freq, paste0(out, ".chr", paste(chrs, collapse="_"), ".freq"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(ci.low, paste0(out, ".chr", paste(chrs, collapse="_"), ".lowCI.freq"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(ci.up, paste0(out, ".chr", paste(chrs, collapse="_"), ".highCI.freq"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
