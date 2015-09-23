## Set up the populations, include.reads, include counts file paths etc, required for
## genome_wide_scan_reads, and other analyses.
## Reads the first two or three arguments. 

########################################################################
## Details
chr <- 1                                #set manually, or from --args
version <- "vx" #v6, v7 etc...
results.tag <- ""

cA <- commandArgs(TRUE)
if(length(cA)){
  chr <- cA[1]
  version <- cA[2]
  if(length(cA)>2){
    results.tag <- cA[3]
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
read.root <- paste0("~/data/", version, "/reads/jj2")
indfile <- paste0("~/data/", version, "/use/", version,"1kg_europe2names.ind")
error.prob <- 0.001

pops <- monocheck <- A <- NA            #To fill in
include.reads <- include.counts <- list()
########################################################################
## Setup according to version.

if(version=="v6" | version=="v7"){

  anc.pops <- c("WHG", "EN", "Yamnaya")
  mod.pops <- c("CEU", "GBR", "IBS", "TSI")
  pops <- c(anc.pops, mod.pops)
#Check if the SNP is monomorphic in these populations. 
  monocheck <- c("CEU", "GBR", "IBS", "TSI", "HungaryGamba_HG", "Loschbour", "Stuttgart",
               "LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN", "Yamnaya")
  A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.227, 0, 0.712, 0.288),3, 4) 

  include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "EN"="Stuttgart",
    "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
  
include.reads <- list(                  #Include these populations as reads
    ## "WHG"=c("LaBrana1", "HungaryGamba_HG"), #Replace LaBrana1 with SpanishMesolithic for the high coverage LaBrana I0585
    "WHG"=c("SpanishMesolithic", "HungaryGamba_HG"), #Replace LaBrana1 with SpanishMesolithic for the high coverage LaBrana I0585
    "EN"=c("LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN"), 
    "Yamnaya"="Yamnaya")
}
if(version=="v7"){
  include.reads[["WHG"]] <- gsub("SpanishMesolithic", "Iberian_Mesolithic", include.reads[["WHG"]], fixed=TRUE)
}
if(version=="v8"){
  mix.dir <- "~/selection/code/files/v8/mixtures/"
  
  if(results.tag==""){stop("Must specify results tag for v8 analysis")}
  include.counts <- list( "CEU"="CEU", "GBR"="GBR", "IBS"="IBS", "TSI"="TSI" )
  always.counts <- c("Loschbour", "Stuttgart")
  
  group <- results.tag
  choice <- read.table(paste0(mix.dir, "Choice", results.tag), as.is=TRUE, header=FALSE)
  include.reads <- list(c(), c(), c())
  names(include.reads) <- unique(choice[,2])
  for(i in 1:NROW(choice)){
    if(choice[i,1] %in% always.counts){
      include.counts[[choice[i,2]]] <- choice[i,1]
    } else{
      include.reads[[choice[i,2]]] <- c(include.reads[[choice[i,2]]], choice[i,1])
    }
  }
  mix.mat <- read.table(paste0(mix.dir, "Proportion", results.tag), as.is=TRUE, header=TRUE)
  rownames(mix.mat) <- mix.mat[,1]
  mix.mat <- mix.mat[,2:NCOL(mix.mat)]
  
  anc.pops <- names(include.reads)
  mod.pops <- c("CEU", "GBR", "IBS", "TSI")
  pops <- c(anc.pops, mod.pops)
  A <- t(mix.mat)[anc.pops,mod.pops]

  ## monocheck <- c(unlist(include.reads), unlist(include.counts))
  ## names(monocheck) <- NULL

  monocheck <-  c("CEU", "GBR", "IBS", "TSI", "YRI")
  
  cat("Monocheck\n")
  print(monocheck)
}
if(version=="peru"){
  mix.dir <- "~/selection/code/files/peru/mixtures/"
  
  if(results.tag==""){stop("Must specify results tag for v8 analysis")}

  
  include.counts <- list( "IBS"="IBS", "YRI"="YRI", "MXL"="MXL", "CLM"="CLM", "PEL"="PEL", "PUR"="PUR")

  always.counts <- c()

  group <- results.tag
  choice <- read.table(paste0(mix.dir, "Choice", results.tag), as.is=TRUE, header=FALSE)
  include.reads <- list(c(), c(), c())
  names(include.reads) <- unique(choice[,2])
  for(i in 1:NROW(choice)){
    if(choice[i,1] %in% always.counts){
      include.counts[[choice[i,2]]] <- choice[i,1]
    } else{
      include.reads[[choice[i,2]]] <- c(include.reads[[choice[i,2]]], choice[i,1])
    }
  }
  mix.mat <- read.table(paste0(mix.dir, "Proportion", results.tag), as.is=TRUE, header=TRUE)
  rownames(mix.mat) <- mix.mat[,1]
  mix.mat <- mix.mat[,2:NCOL(mix.mat)]

  anc.pops <- colnames(mix.mat)
  mod.pops <- rownames(mix.mat)
  
  pops <- c(anc.pops, mod.pops)
  A <- t(mix.mat)[anc.pops,mod.pops,drop=FALSE]

  ## monocheck <- c(unlist(include.reads), unlist(include.counts))
  ## names(monocheck) <- NULL

  monocheck <-  names(include.counts)
  
  cat("Monocheck\n")
  print(monocheck)
}

