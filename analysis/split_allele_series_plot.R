## Plot the time series of an alelle split between the polygenic
## selection populations. 

library(plotrix)
library(RColorBrewer)

########################################################################
## Details
root <- "~/selection/counts/all"
readmefile <- "~/selection/analysis/series/figure2.readme"
out <- "~/selection/analysis/series/"
ang <- 20
ylim <- c(0,1)

########################################################################
## POPS

pp <- read.table("~/selection/analysis/series/polypop_list.txt", as.is=TRUE)

########################################################################
## Details
snpname <- "rs12913832"
flip <- TRUE
outname <- "herc2series.pdf"

data <- read.table("~/selection/counts/all.reads.freq", as.is=TRUE, header=TRUE)
freq <- data[,6:NCOL(data)]
data <- data[,1:5]
rownames(freq) <- data[,1]

pops <- unique(pp[,1])

ff<-freq[snpname,]
if(flip){ff <- 1-ff}
pp$f <- unlist(ff[pp[,2]])
         
cols <- brewer.pal( length(pops), "Set1")

pdf(paste0(out, outname), width=12, height=6)
plot(0,0, col="white", xlim=c(-8000, 0), ylim=ylim, bty="n", xlab="Years before present", ylab="HERC2 allele frequency")

for(i in 1:length(pops)){
    pop <- pops[i]
    
    int.starts <- pp[pp[,1]==pop,3]
    int.ends <- pp[pp[,1]==pop,4]
    af <- pp[pp[,1]==pop,5]

    rect(-int.starts, af-0.02, -int.ends, af+0.02, col=cols[i], density=20, angle=ang*i, border=cols[i])
    ln <- length(int.ends)

    if(ln>1){
        segments(-int.ends[1:(ln-1)], af[1:(ln-1)], -int.starts[2:ln], af[2:ln], col=cols[i], lwd=2)
    }
}
legend("bottomright", pops, col=cols,lwd=2, bty="n")

dev.off()
