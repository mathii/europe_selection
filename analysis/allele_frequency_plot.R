## box series plots.
## Plot the series for individual snps. 

source("~/selection/code/lib/readlib.R")
source("~/selection/code/lib/plotlib.R")
  
library(plotrix)
library(RColorBrewer)

########################################################################
## Details
root <- "~/selection/counts/all"
readmefile <- "~/selection/analysis/series/figure2.readme"
out <- "~/selection/analysis/series/"
ang <- 12
ylim <- c(0,1)
outsize=c(6,6)



## ########################################################################

## int.names <- c("IBS", "TSI", "CEU", "GBR", "FIN")
## long.names <- c("Spanish", "Italian", "Central", "British", "Finnish")
## int.starts <- c(10000, 8000, 6000, 4000, 2000)
## int.ends <- c(9000, 7000, 5000, 3000, 1000)
## int.include <- c("IBS", "TSI", "CEU", "GBR", "FIN")
## names(int.include) <-c("IBS", "TSI", "CEU", "GBR", "FIN")
## outname <- "figure2c.pdf"

## ########################################################################

int.names <- c("SHG","Yamnaya")
long.names <-  c("Swedish Hunter Gatherers",  "Yamnaya")
int.starts <- c(10000, 4000)
int.ends <- c(6000, 0)
int.include <-  c("SHG", "SHG", "SHG", "Yamnaya")
names(int.include) <-  c("SwedenSkoglund_MHG", "Motala_HG", "SwedenSkoglund_NHG", "Yamnaya")
outname <- "figure2b.pdf"

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

pdf(paste0(out, outname), width=6, height=6)
plot(0,0, col="white", xlim=c(-max(int.starts), 0), ylim=ylim, bty="n", xlab="", ylab="Frequency", xaxt="n")

for(i in 1:NROW(data)){
    int.counts <- int.totals <- rep(0, length(int.names))
    names(int.counts) <- names(int.totals) <- int.names
    for(j in 1:NCOL(counts)){
        if(!(names(counts)[j] %in% names(int.include))){
            cat(paste0("Population ", names(counts[j]), " not included\n"))
            next
        }
        int.counts[int.include[names(counts)[j]]] <- int.counts[int.include[names(counts)[j]]]+counts[i,j]
        int.totals[int.include[names(totals)[j]]] <- int.totals[int.include[names(totals)[j]]]+totals[i,j]
    }

    
    f=int.counts/int.totals
    if(readme[data[i,1],3]==data[i,4]){f=1-f}
    rect(-int.starts, f-0.01, -int.ends, f+0.01, col=cols[i], density=20, angle=ang*i, border=cols[i])
    ln <- length(int.ends)

    text(-0.5*(int.ends+int.starts)+100*i-50*ln, f+ifelse(f<0.7,0.03,-0.03), as.character(int.totals), col=cols[i], cex=0.7)
    for(j in 1:NCOL(counts)){
        
    }
}

mtext(long.names, side=3, at=-0.5*(int.ends+int.starts), cex=0.8)
dev.off()
