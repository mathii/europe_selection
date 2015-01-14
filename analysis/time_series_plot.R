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
outname <- "figure2a.pdf"
what <- NA
ang <- 20
ylim <- c(0,1)

########################################################################

int.names <- c("WHG", "EN", "MN", "LN/BA","CEU")
long.names <- c("Western Hunter\nGatherers", "Early neolithic", "Middle neolithic", "Late neolithic/Bronze age", "CEU")
int.starts <- c(8000, 7200, 5800, 4800, 100)
int.ends <- c(7700, 6900, 5200, 3600, -50)
int.include <- c("WHG", "WHG", "WHG", "EN", "EN", "EN", "EN", "EN", "EN", "MN", "MN", "MN", "MN", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "CEU")
names(int.include) <-c("Loschbour", "LaBrana1", "HungaryGamba_HG", "Starcevo_EN", "Stuttgart", "Spain_EN", "LBK_EN", "LBKT_EN", "HungaryGamba_EN", "Spain_MN", "Baalberge_MN", "Iceman", "Esperstedt_MN", "HungaryGamba_CA", "Alberstedt_LN", "Corded_Ware_LN", "Bell_Beaker_LN", "BenzigerodeHeimburg_LN", "Unetice_EBA", "HungaryGamba_BA", "Halberstadt_LBA", "CEU")

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
    cols <- c(brewer.pal( 9, "Set1"), rainbow(NROW(data)-9))
} else{
    cols <- brewer.pal( NROW(data), "Set1")
    cols[6] <- "darkgrey"
}

pdf(paste0(out, outname), width=12, height=6)
plot(0,0, col="white", xlim=c(-max(int.starts), 0), ylim=ylim, bty="n", xlab="Years", ylab="Frequency")

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

    include.segments <- rep(TRUE, ln-1)
    for(k in 1:(ln-1)){
        mat <- matrix(c(int.counts[k], int.totals[k]-int.counts[k], int.counts[k+1], int.totals[k+1]-int.counts[k+1]),2,2)
        include.segments[k] <- fisher.test(mat)$p.val<0.05
    }
    
    segments(-int.ends[1:(ln-1)][include.segments], f[1:(ln-1)][include.segments], -int.starts[2:ln][include.segments], f[2:ln][include.segments], col=cols[i], lwd=2)
    segments(-int.ends[1:(ln-1)][!include.segments], f[1:(ln-1)][!include.segments], -int.starts[2:ln][!include.segments], f[2:ln][!include.segments], col=cols[i], lwd=2, lty=2)

    text(-0.5*(int.ends[1:(ln-1)]+int.starts[1:(ln-1)])+100*i-50*ln, f[1:(ln-1)]+ifelse(f[1:(ln-1)]<0.7,0.03,-0.03), as.character(int.totals[1:(ln-1)]), col=cols[i], cex=0.7)
    for(j in 1:NCOL(counts)){
        
    }
}
## legend("bottomright", paste0(readme[,1], " (0=", readme[,4], ", 1=", readme[,6], ")"), col=cols,lwd=2, bty="n", cex=0.75)
legend("bottomright", paste0(readme[,1], " (", readme[,2], ")"), col=cols,lwd=2, bty="n", cex=0.75)

mtext(long.names, side=3, at=-0.5*(int.ends+int.starts), cex=0.8)
dev.off()
