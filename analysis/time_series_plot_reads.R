## box series plots.
## Plot the series for individual snps. 

source("~/selection/code/lib/readlib.R")
source("~/selection/code/lib/plotlib.R")
source("~/selection/code/lib/3pop_lib.R")

library(plotrix)
library(RColorBrewer)

########################################################################
## Details
outname <- "figure2a.pdf"
what <- NA
version <- NA
ang <- 20
ylim <- c(0,1)
error.prob <- 0.001
readmefile <- "~/selection/code/files/figure2.readme" 

########################################################################

if(length(commandArgs(TRUE))){
    what <- commandArgs(TRUE)[1]
    version <- commandArgs(TRUE)[2]
    readmefile <- paste0("~/selection/code/files/", what, ".readme")
    outname <- paste0(what, "a.pdf")
}

root <- paste0("~/selection/counts/", version, "/all")
read.root <- paste0("~/data/",version,"/reads/")
indfile <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.ind")
out <- paste0("~/selection/analysis/",version,"/series/")


########################################################################

int.names <- c("WHG", "EN", "MN", "LN/BA","CEU/IBS")
long.names <- c("Western Hunter\nGatherers", "Early neolithic", "Middle neolithic", "Late neolithic/Bronze age", "CEU/IBS")
int.starts <- c(8000, 7200, 5800, 4800, 100)
int.ends <- c(7700, 6900, 5200, 3600, -50)
int.include <- c("WHG", "WHG", "WHG", "EN", "EN", "EN", "EN", "EN", "EN", "MN", "MN", "MN", "MN", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "CEU/IBS", "CEU/IBS")
names(int.include) <-c("Loschbour", "Iberian_Mesolithic", "HungaryGamba_HG", "Starcevo_EN", "Stuttgart", "Spain_EN", "LBK_EN", "LBKT_EN", "HungaryGamba_EN", "Spain_MN", "Baalberge_MN", "Iceman", "Esperstedt_MN", "HungaryGamba_CA", "Alberstedt_LN", "Corded_Ware_LN", "Bell_Beaker_LN", "BenzigerodeHeimburg_LN", "Unetice_EBA", "HungaryGamba_BA", "Halberstadt_LBA", "CEU", "IBS")

include.reads <- list(                  #Include these populations as reads
    "WHG"=c("Iberian_Mesolithic", "HungaryGamba_HG"), #SpanishMesolithic is the high coverage LaBrana1 I0585 and LaBrana2
    "EN"=c("LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN"),
    "MN"=c( "Spain_MN", "Baalberge_MN", "Iceman", "Esperstedt_MN" ),
    "LN/BA"=c( "HungaryGamba_CA", "Alberstedt_LN", "Corded_Ware_LN", "Bell_Beaker_LN", "BenzigerodeHeimburg_LN", "Unetice_EBA", "HungaryGamba_BA", "Halberstadt_LBA")
    )
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "EN"="Stuttgart",
    "CEU/IBS"=c("CEU", "IBS")
     )

## ########################################################################

## Setup the count data
counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- counts[,6:NCOL(counts)]
totals <- totals[,6:NCOL(totals)]

## Readme data and filter count data
readme <- read.table(readmefile, as.is=TRUE)
rownames(readme) <- readme[,2]
rownames(data) <- data[,1]
rownames(counts) <- data[,1]
rownames(totals) <- data[,1]
data <- data[rownames(readme),]
counts <- counts[rownames(readme),]
totals <- totals[rownames(readme),]

if(NROW(data)>9){
    cols <- c(brewer.pal( 9, "Set1"), rainbow(NROW(data)-9))
} else{
    cols <- brewer.pal( NROW(data), "Set1")
    cols[6] <- "darkgrey"
}

empty.data <- make.empty.data(int.names)
include.read.samples <- read.samples(indfile, include.reads)
    
pdf(paste0(out, outname), width=12, height=6)
plot(0,0, col="white", xlim=c(-max(int.starts), 0), ylim=ylim, bty="n", xlab="Years", ylab="Frequency")

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

    include.segments <- rep(TRUE, ln-1)
    ## include.segments <- c(FALSE, TRUE, FALSE, TRUE)
    ## TODO - include segments based on p value? 
    
    segments(-int.ends[1:(ln-1)][include.segments], f[1:(ln-1)][include.segments], -int.starts[2:ln][include.segments], f[2:ln][include.segments], col=cols[i], lwd=2)
    ## segments(-int.ends[1:(ln-1)][!include.segments], f[1:(ln-1)][!include.segments], -int.starts[2:ln][!include.segments], f[2:ln][!include.segments], col=cols[i], lwd=2, lty=2)

    ## text(-0.5*(int.ends[1:(ln-1)]+int.starts[1:(ln-1)])+100*i-50*ln, f[1:(ln-1)]+ifelse(f[1:(ln-1)]<0.7,0.03,-0.03), format(eff.totals[1:(ln-1)], digits=2, nsmall=1), col=cols[i], cex=0.7)
    text(-0.5*(int.ends[1:(ln-1)]+int.starts[1:(ln-1)]), f[1:(ln-1)]+ifelse(f[1:(ln-1)]<0.7,0.03,-0.03), format(eff.totals[1:(ln-1)], digits=2, nsmall=1), col=cols[i], cex=0.7)

}
## legend("bottomright", paste0(readme[,1], " (0=", readme[,4], ", 1=", readme[,6], ")"), col=cols,lwd=2, bty="n", cex=0.75)
legend("bottomright", paste0(readme[,1], " (", readme[,2], ")"), col=cols,lwd=2, bty="n", cex=0.75)

mtext(long.names, side=3, at=-0.5*(int.ends+int.starts), cex=0.8)
dev.off()
