## Estimate the distribution of genetic height, for the ancient samples
library(truncnorm)

## Use the read count information to estimate the freqeuncies in each population
## Should really be implemented in spindrift. Actually this just counts
## Here we're using the polymap populations. 
source("~/selection/code/lib/3pop_lib.R")
source("~/selection/code/lib/readlib.R")

########################################################################
## Details
verbose=TRUE
which.map <- ""
if( Sys.info()["login"]!=Sys.info()["user"]){
    verbose=FALSE
}

MIN.N.CHR <- 3

########################################################################
#Arguments

cA <- commandArgs(TRUE)
if(!length(cA)!=2){
    stop("Must specify version and snplist")
}
version <- cA[1]
snplist <- cA[2]
pops.to.use <- NA
extag <- ""
if(length(cA)>3){
    pops.to.use <- strsplit(cA[3], ",", fixed=TRUE)[[1]]
    extag <- paste0("_",gsub(",", "_", cA[3], fixed=TRUE))
}
########################################################################
#MCMC parameters
burn.in <- 100
N.iter <- 1000
if(length(cA)>2){
    burn.in <- as.numeric(cA[4])
    N.iter <- as.numeric(cA[5])
}
thin <- 1                               #Why bother?

if(length(cA)>5
){
    MIN.N.CHR <- cA[6]
}

if(length(cA)>6
){
    which.map <- cA[7]
}

########################################################################
## Details
root <- paste0("~/selection/counts/",version,"/all")
read.root <- paste0("~/data/",version,"/reads/jj2")
indfile <- paste0("~/data/",version,"/use/", version, "1kg_europe2names.ind")
snpfile <- paste0("~/data/",version,"/use/",version ,"1kg_europe2names.snp")
polymap <- paste0("~/selection/code/files/",version,"/polymap.txt")
gwas <- read.table(paste0("~/selection/data/gwas/",snplist,".gwas"), as.is=TRUE)
colnames(gwas) <- c("CHR", "POS", "EFFECT", "OTHER", "BETA")
if(which.map!=""){
    polymap <- paste0("~/selection/code/files/",version ,"/polymap.", which.map, ".txt" )
}

error.prob <- 0.01
########################################################################
#Include these populations as hard calls
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "Central_EN"="Stuttgart",
    "CEU"="CEU",  "IBS"="IBS" )
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
    "HG"="Loschbour",
    "CEM"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
}

if(version=="sard"){
  include.counts <- list(                 #Include these populations as hard calls. 
    "HG"="Loschbour",
    "CEM"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI", "Sardinian"="Sardinian" )
}

if(version=="bell"){
include.counts <- list(                 #Include these populations as hard calls. 
    "CEU"="CEU", "IBS"="IBS" )
}


########################################################################
#Standard setup

include.reads <- list()
polymap <- read.table(polymap, as.is=TRUE, header=FALSE)
for(i in 1:NROW(polymap)){
    if(polymap[i,2] %in% exclude){next}
    
    if(polymap[i,2] %in% names(include.reads)){
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

data$CHR.POS <- paste0(data$CHR, "_", data$POS)
gwas$CHR.POS <- paste0(gwas$CHR, "_", gwas$POS)
include <- data$CHR.POS %in% gwas$CHR.POS
counts <- counts[include,]
totals <- totals[include,]
data <- data[include,]
gwas.include <- gwas$CHR.POS %in% data$CHR.POS
gwas <- gwas[gwas.include,]
if(!all(data$CHR.POS==gwas$CHR.POS)){stop("Data and gwas not matched")}

## get list of samples in each population of reads
include.read.samples <- read.samples(indfile, include.reads, c(exclude, unlist(include.counts)))

########################################################################
#Flip gwas alleles
#So that Betas are for the REFERENCE allele relative to the ALT allele
for(i in 1:NROW(gwas)){
    if(all(gwas[i,c("OTHER", "EFFECT")]==data[i,c("REF", "ALT")])){ #If the ALT allele is the effect allele
        gwas[i,"BETA"] <- -gwas[i,"BETA"]
        gwas[i,c("OTHER", "EFFECT")] <- gwas[i,c("EFFECT", "OTHER")] 
    }

}
    
########################################################################
#Get population frequencies for each snp

empty.data <- make.empty.data(pops)
freq.data.list <- list()

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

    
    this.snp <- data[i,1]
    ## this.read <- reads[reads[,1]==this.snp,]
    this.read <- reads[(1+N.read.samples*(readi-1)):(N.read.samples*readi),]
    if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

    freq.data.list[[i]] <- make.freq.data(pops, include.reads, include.read.samples, include.counts,
                                this.read, counts[i,], totals[i,], empty.data)

    readi <- readi+1
}

########################################################################
## SNPs with sufficient data
data.size <-  data.frame(do.call(rbind, lapply(freq.data.list, effective.data.size)))
include <- apply(data.size, 1, min)>MIN.N.CHR
freq.data.list <- freq.data.list[include]
data <- data[include,]
gwas <- gwas[include,]

########################################################################
#Mean genetic value estimate
freq.est <- matrix(NA, nrow=NROW(data), ncol=length(pops))
colnames(freq.est) <- pops
for(i in 1:NROW(data)){
  fr <- fit.unconstrained.model.reads(freq.data.list[[i]], error.prob=error.prob)$par
  ## Replace missing values with the mean of the non-missing values. Shrinking towards mean
  fr[is.na(fr)] <- mean(fr, na.rm=TRUE)
  freq.est[i,] <- fr
}

fr.4 <- round(colSums(freq.est*gwas$BETA),4)

## dsrs <- read.table("~/selection/analysis/v6/poly/Height.old/pred_height.txt", as.is=TRUE)
## dsrs[,5] <- fr.4[dsrs[,2]]
## write.table(dsrs, "~/selection/analysis/v6/poly/Height/pred_height.txt", col.names=FALSE, row.names=FALSE, sep ="\t", quote=FALSE)

########################################################################
#MCMC posterior interval estimate.

## proposal <- function(x, sd=0.001){return(rtruncnorm(length(x), mean=x, sd=sd, a=0, b=1))}
## proposal <- function(x, sd=0.001){return(rnorm(length(x), mean=x, sd=sd))}
proposal <- function(x, sd=0.001){
    prop <- rnorm(length(x), mean=x, sd=sd)
    resample <- (prop<0 | prop>1)
    if(any(resample)){
        prop[resample] <- rtruncnorm(sum(resample), mean=x[resample], sd=sd, a=0, b=1)
    }
    return(prop)
}

## Implicit uniform (0,1) prior on frequencies
mh.step <- function(f, prop.freq.data.list, proposal.sd=0.001){
    new.f <- proposal(f, sd=proposal.sd)

    if(any(new.f>1)|any(new.f<0)){
        ratio <- 0
    }else{
        old.ll <- sum(mapply(likelihood.reads, f, pop.freq.data.list))
        new.ll <- sum(mapply(likelihood.reads, new.f, pop.freq.data.list))
        #Normalising constant because we are truncating the proposal density. 
        trunc.const <- sum(log((pnorm(1,mean=f,sd=proposal.sd)-pnorm(0,mean=f,sd=proposal.sd))/(pnorm(1,mean=new.f,sd=proposal.sd)-pnorm(0,mean=new.f,sd=proposal.sd))))
        ratio <- exp(new.ll-old.ll+trunc.const)
    }
    if(ratio>1 | runif(1)<ratio){
        return(list(f=new.f, accept=TRUE))
    }else{
        return(list(f=f, accept=FALSE))
    }
}

if(all(is.na(pops.to.use))){pops.to.use <- pops}

results <- matrix(0, nrow=length(pops.to.use), ncol=4)
rownames(results) <- pops.to.use
colnames(results) <- c("MLE", "Post.Mn", "Post.5", "Post.95")
results[,1] <- fr.4[pops.to.use]

popi <- 1
for(pop in pops.to.use){
    prop.sd <- 0.01
    if(pop %in% c("CEU", "IBS")){prop.sd <- 0.003}
    
    cat(paste0(pop, "\n"))
    pop.freq.data.list <- list()
    for(i in 1:NROW(data)){
        pop.freq.data.list[[i]] <- list(freq.data.list[[i]][[pop]])
        names(pop.freq.data.list[[i]]) <- pop
    }

    #Eek
    ## f <- unlist(sapply(pop.freq.data.list, fit.unconstrained.model.reads, error.prob=error.prob)[1,])
    f <- runif(NROW(data))
    
    cat("Burn in...")
    accept <- c(0,0)
    for(i in 1:burn.in){
        if(verbose){cat(paste0("\rBurn in...", i))}
        new.res <- mh.step(f, pop.freq.data.list, proposal.sd=prop.sd)
        f <- new.res$f
        accept <- accept + c(new.res$accept, 1)
    }
    cat(paste0("\rBurn in...Done; Acceptance ratio: ", accept[1]/accept[2], "\n"))

    cat("Starting iterations...")
    accept <- c(0,0)
    h.values <- rep(NA, N.iter/thin)
    for(i in 1:N.iter){
        if(verbose){cat(paste0("\rStarting iterations...", i))}
        new.res <- mh.step(f, pop.freq.data.list, proposal.sd=prop.sd)
        f <- new.res$f
        accept <- accept + c(new.res$accept, 1)
        if(!(i%%thin)){
            h.values[i/thin] <- sum(f*gwas$BETA)
        }
    }
    cat(paste0("\rStarting iterations...Done; Acceptance ratio: ", accept[1]/accept[2], "\n"))

    results[popi,2:4] <- c(mean(h.values), as.numeric(quantile(h.values, c(0.05, 0.95))))
    popi <- popi+1
}

write.table(results, paste0("~/selection/analysis/", version, "/series/height_series_mcmc_estimates", extag,".",snplist, ".minchr", MIN.N.CHR, ifelse(which.map=="", "", paste0(".", which.map)), ".txt"), col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t") 



