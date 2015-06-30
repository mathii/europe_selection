#This is the genome-wide scan but using the read level information to try and
#get some idea about diploid calls. For reasons of speed, we break this one
#by chromosome. Genomic correction has to be done in another step for this reason.

#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
#Using genotype probabilities imputed with Beagle.

source("~/selection/code/lib/3pop_lib.R")

#Modern GBR, CEU, IBS, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
chr <- 1                                #set manually, or from --args
version <- "vx" #v6, v7 etc...
results.tag <- ""
which.impute <- "within"

cA <- commandArgs(TRUE)
if(length(cA)){
  chr <- cA[1]
  version <- cA[2]
  if(length(cA)>2){
    results.tag <- cA[3]
  }
  if(length(cA)>3){
    which.impute <- cA[4]
  }
}

verbose=TRUE
## Supposed to check if running on cluster, but YMMV
if( Sys.info()["login"]!=Sys.info()["user"]){
    verbose=FALSE
}

########################################################################
## Details
root <- paste0("~/selection/counts/", version, "/all")
out <- paste0("~/selection/analysis/", version, "/gscan/")
indfile <- paste0("~/data/", version, "/use/", version,"1kg_europe2names.ind")
impute.file <- paste0("~/selection/imputation/", version, "/imputed.", which.impute ,".chr", chr, ".vcf.gz")
error.prob <- 0.001

pops <- c("WHG", "EN", "Yamnaya", "CEU", "GBR", "IBS", "TSI")
#Check if the SNP is monomorphic in these populations. 
monocheck <- c("CEU", "GBR", "IBS", "TSI", "HungaryGamba_HG", "Loschbour", "Stuttgart",
               "LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN", "Yamnaya")
A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.227, 0, 0.712, 0.288),3, 4) 

########################################################################

if(version=="v6" | version=="v7"){

  pops <- c("WHG", "EN", "Yamnaya", "CEU", "GBR", "IBS", "TSI")
#Check if the SNP is monomorphic in these populations. 
  monocheck <- c("CEU", "GBR", "IBS", "TSI", "HungaryGamba_HG", "Loschbour", "Stuttgart",
               "LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN", "Yamnaya")
  A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.227, 0, 0.712, 0.288),3, 4) 

  include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "EN"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
  
include.probs <- list(                  #Include these populations as reads
    ## "WHG"=c("LaBrana1", "HungaryGamba_HG"), #Replace LaBrana1 with SpanishMesolithic for the high coverage LaBrana I0585
    "WHG"=c("SpanishMesolithic", "HungaryGamba_HG"), #Replace LaBrana1 with SpanishMesolithic for the high coverage LaBrana I0585
    "EN"=c("LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN"), 
    "Yamnaya"="Yamnaya")
}
if(version=="v7"){
  include.probs[["WHG"]] <- gsub("SpanishMesolithic", "Iberian_Mesolithic", include.probs[["WHG"]], fixed=TRUE)
}
if(version=="v8"){
  mix.dir <- "~/selection/code/files/v8/mixtures/"
  
  if(results.tag==""){stop("Must specify results tag - group from 1-6 - for v8 analysis")}
  include.counts <- list( "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
  always.counts <- c("Loschbour", "Stuttgart")
  group <- results.tag
  choice <- read.table(paste0(mix.dir, "Choice", results.tag), as.is=TRUE, header=FALSE)
  include.probs <- list(c(), c(), c())
  names(include.probs) <- unique(choice[,2])
  for(i in 1:NROW(choice)){
    if(choice[i,1] %in% c("Loschbour", "Stuttgart")){
      include.counts[[choice[i,2]]] <- choice[i,1]
    } else{
      include.probs[[choice[i,2]]] <- c(include.probs[[choice[i,2]]], choice[i,1])
    }
  }
  mix.mat <- read.table(paste0(mix.dir, "Proportion", results.tag), as.is=TRUE, header=TRUE)
  rownames(mix.mat) <- mix.mat[,1]
  mix.mat <- mix.mat[,2:NCOL(mix.mat)]
  
  anc.pops <- names(include.probs)
  mod.pops <- c("CEU", "GBR", "IBS", "TSI")
  pops <- c(anc.pops, mod.pops)
  A <- t(mix.mat)[anc.pops,mod.pops]

  monocheck <- c(unlist(include.probs), unlist(include.counts))
  names(monocheck) <- NULL
}

##################################################################################################

## Setup the count data. 
counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
include <- data$CHR==chr
data <- data[include,]

counts <- data.matrix(counts[,6:NCOL(counts)])
totals <- data.matrix(totals[,6:NCOL(totals)])
counts <- counts[include,]
totals <- totals[include,]

#Load imputed likelihoods
impute <- read.table(impute.file, comment.char="", as.is=TRUE, header=FALSE, sep="\t", fill=TRUE)
comment.lines <- sum(grepl("^##", impute[,1]))
impute <- read.table(impute.file, comment.char="", as.is=TRUE, header=TRUE, skip=comment.lines, sep="\t", fill=TRUE)
rownames(impute) <- impute$ID
impute <- impute[impute$ID %in% data[,1],]
impute.info <- impute[,8]
impute <- impute[,10:NCOL(impute)]

## get list of samples in each population of reads
include.prob.samples <- read.samples(indfile, include.probs)

## set up results
results <- matrix(0, nrow=NROW(data), ncol=3)
rownames(results) <- data$ID

## Data structure
empty.data <- make.empty.data(pops)

for(i in 1:NROW(data)){
    this.snp <- data[i,1]
    if(verbose){cat(paste0("\r", i, " ", this.snp))}
    this.prob <- impute[i,]

    freq.data <- make.prob.freq.data(pops, include.probs, include.prob.samples, include.counts,
                                this.prob, counts[i,], totals[i,], empty.data)
    monomorphic <- all(counts[i,monocheck]==0)|all(counts[i,monocheck]==totals[i,monocheck])
    if(monomorphic){
        results[i,] <- NA
    }else{
        AR2 <- as.numeric(strsplit(strsplit(impute.info[i], ";", fixed=TRUE)[[1]][1], "=", fixed=TRUE)[[1]][2])
        results[i,] <- c(test.3pop.reads(freq.data, A, error.prob=error.prob), AR2)
    }
}

results <- results[!is.na(results[,2]),]

results <- cbind(rownames(results), results)
colnames(results) <- c("ID", "ChiSq", "uncorrected.p", "AR2")
results <- data.frame(results)
out.file <-  paste0("~/selection/analysis/",version,"/gscan/scan_results_imputed", results.tag, ".", which.impute, ".chr", chr, ".txt")
print(out.file)
write.table(results,out.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

