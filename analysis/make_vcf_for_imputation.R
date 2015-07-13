## Does what it says...

source("~/selection/code/lib/3pop_lib.R")
source("~/selection/code/lib/readlib.R")

########################################################################
## Details
chr <- NA                                #set manually, or from --args
root <- ""
cA <- commandArgs(TRUE)
if(length(cA)){
    chr <- as.numeric(cA[1])
    version <- cA[2]
}
data.root <- paste0("~/data/", version, "/use/", version, "1kg_europe2names")

verbose=TRUE
## Supposed to check if running on cluster, but YMMV
if( Sys.info()["login"]!=Sys.info()["user"]){
    verbose=FALSE
}

########################################################################
## Details
root <- paste0("~/selection/counts/",version,"/all")
read.root <- paste0("~/data/",version,"/reads/jj2")
indfile <- paste0(data.root, ".ind")
snpfile <- paste0(data.root, ".snp")
out <- paste0("~/selection/imputation/",version,"/pre.chr", chr, ".vcf")
populations <- scan(paste0("~/selection/code/files/",version,"/used_all.txt"), "")
error.prob <- 0.01

########################################################################
#Include these samples as hard calls
include.hard <- c("Loschbour", "Stuttgart")
#Exclude these
exclude <- c("LaBrana1", "CEU", "FIN", "GBR", "IBS", "TSI", "YRI", "Loschbour", "Stuttgart")

########################################################################

include.reads <- list()
ind <- read.table(indfile, header=FALSE, as.is=TRUE)
snp <- read.table(snpfile, header=FALSE, as.is=TRUE)
snp <- snp[snp[,2]==chr,]

samples <- sort(ind[ind[,3]%in%populations & !(ind[,3]%in%exclude),1])
log10e <- log10(exp(1))

cat("##fileformat=VCFv4.0\n", file=out)
cat("##source=make_vcf_for_imputation.R\n", file=out, append=TRUE)
cat("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", file=out, append=TRUE)
cat("##FORMAT=<ID=GL,Number=3,Type=Integer,Description=\"Genotype Likelihood\">\n", file=out, append=TRUE)
cat(paste(c("#CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples, "\n"), collapse="\t"), file=out, append=TRUE)

reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
reads <- reads[order(match(reads[,1],snp[,1])),]
read.sample.counts <- table(reads[,1])  #Number of samples for each SNP
N.read.samples <- as.numeric(read.sample.counts[1])
if(!all(read.sample.counts==N.read.samples)){stop("Different number of read samples for some snps")}
readi <- 1
cat(paste0("Done\n"))

for(i in 1:NROW(snp)){
    this.gtgl.str <- rep("", length(samples))
    ## names(this.gtgl.str) <- samples

    this.snp <- snp[i,1]
    this.read <- reads[(1+N.read.samples*(readi-1)):(N.read.samples*readi),]
    if(!all(this.read[,1]==this.snp)){stop("Selected the wrong SNP")}

    for(j in 1:length(samples)){
        sample=samples[j]
        ref.alt <- this.read[this.read[,2]==sample,3:4]
        N <- as.numeric(ref.alt[1]+ref.alt[2])
        X <- as.numeric(ref.alt[1])
        p00 <- dbinom(X, N, 1-error.prob, log=TRUE)*log10e
        p01 <- dbinom(X, N, 0.5, log=TRUE)*log10e
        p11 <- dbinom(X, N, error.prob, log=TRUE)*log10e
        this.gtgl.str[j] <- paste0("./.:", format(p00, scientific=FALSE, nsmall=4, digits=4), ",", format(p01, scientific=FALSE, nsmall=4, digits=4), ",", format(p11, scientific=FALSE, nsmall=4, digits=4))
        
    }
    cat(paste(c(chr,snp[i,4],snp[i,1], snp[i,5], snp[i,6], 0, "PASS", ".", "GT:GL", this.gtgl.str, "\n"), collapse="\t"), file=out, append=TRUE)

    readi=readi+1
}
