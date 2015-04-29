## box series plots.
## Plot the series for individual snps. 

source("~/selection/code/lib/readlib.R")
source("~/selection/code/lib/plotlib.R")
source("~/selection/code/lib/3pop_lib.R")

library(plotrix)
library(RColorBrewer)

########################################################################
## Details
ang <- 12
ylim <- c(0,1)
outsize=c(6,6)
error.prob <- 0.001

what <- "figure2"
sub="c"
version <- ""
if(length(commandArgs(TRUE))){
    what=commandArgs(TRUE)[1]
    sub=commandArgs(TRUE)[2]
    version <- commandArgs(TRUE)[3]
}

outname <- paste0(what, sub, ".pdf")
readmefile <- paste0("~/selection/code/files/", what, ".readme")
root <- paste0("~/selection/counts/",version,"/all")
read.root <- paste0("~/data/",version,"/reads/")
indfile <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.ind")
out <- paste0("~/selection/analysis/",version,"/series/")

## ########################################################################

int.names <- c("IBS", "TSI", "CEU", "GBR", "FIN")
long.names <- c("Spanish", "Italian", "Central", "British", "Finnish")
int.starts <- c(10000, 8000, 6000, 4000, 2000)
int.ends <- c(9000, 7000, 5000, 3000, 1000)
int.include <- c("IBS", "TSI", "CEU", "GBR", "FIN")
names(int.include) <-c("IBS", "TSI", "CEU", "GBR", "FIN")
include.counts <- list("IBS"="IBS", "TSI"="TSI", "GBR"="GBR", "CEU"="CEU", "FIN"="FIN")
include.reads <- list()

## ########################################################################
if(sub=="b"){
    int.names <- c("SHG","Yamnaya")
    long.names <-  c("Swedish Hunter Gatherers",  "Yamnaya")
    int.starts <- c(10000, 4000)
    int.ends <- c(6000, 0)
    int.include <-  c("SHG", "SHG", "SHG", "Yamnaya")
    names(int.include) <-  c("SwedenSkoglund_MHG", "Motala_HG", "SwedenSkoglund_NHG", "Yamnaya")
    include.reads <- list( "SHG"=c("SwedenSkoglund_MHG", "Motala_HG", "SwedenSkoglund_NHG"),
                          "Yamnaya"="Yamnaya")
    include.counts <- list()
}

## ########################################################################

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- counts[,6:NCOL(counts)]
totals <- totals[,6:NCOL(totals)]

readme <- read.table(readmefile, as.is=TRUE)
rownames(readme) <- readme[,2]
rownames(data) <- data[,1]
rownames(counts) <- data[,1]
rownames(totals) <- data[,1]

data <- data[rownames(readme),]
counts <- counts[rownames(readme),]
totals <- totals[rownames(readme),]

if(NROW(data)>9){
    cols <- c(brewer.pal( 9, "Set2"), rainbow(NROW(data)-9))
} else{
    cols <- brewer.pal( NROW(data), "Set1")
    cols[6] <- "darkgrey"
}

empty.data <- make.empty.data(int.names)
include.read.samples <- read.samples(indfile, include.reads)

pdf(paste0(out, outname), width=6, height=6)
plot(0,0, col="white", xlim=c(-max(int.starts), 0), ylim=ylim, bty="n", xlab="", ylab="Frequency", xaxt="n")

for(i in 1:NROW(data)){
    cat(paste0("\r", readme[i,2]))
    chr <- data[i,"CHR"]
    read.data <- read.table(paste0(read.root, "jj2.chr", chr, ".readcounts.gz"), as.is=TRUE)
    read.data <- read.data[read.data[,1]==data[i,"ID"],]
    
    freq.data <- make.freq.data(int.names, include.reads, include.read.samples, include.counts,
                                read.data, counts[i,], totals[i,], empty.data)
    f <- fit.unconstrained.model.reads(freq.data, error.prob=error.prob)$p
    names(f) <- int.names
    eff.totals <- round(effective.data.size(freq.data),1)
    
    if(readme[data[i,1],3]==data[i,4]){f=1-f}
    rect(-int.starts, f-0.01, -int.ends, f+0.01, col=cols[i], density=20, angle=ang*i, border=cols[i])
    ln <- length(int.ends)

    text(-0.5*(int.ends+int.starts), f+ifelse(f<0.7,0.03,-0.03), format(eff.totals, digits=2, nsmall=1), col=cols[i], cex=0.7)
    for(j in 1:NCOL(counts)){
        
    }
}

mtext(long.names, side=3, at=-0.5*(int.ends+int.starts), cex=0.8)
dev.off()
