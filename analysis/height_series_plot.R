## box series plot of height. 
## Just going to hardcode the results in here for now.

library(plotrix)
library(RColorBrewer)

########################################################################
## Details
root <- "~/selection/counts/all"
readmefile <- "~/selection/analysis/series/figure2.readme"
out <- "~/selection/analysis/poly/Height/"
outname <- "heighseries.pdf"
ang <- 20
ylim <- c(-1,1)

########################################################################
## Details
data <- read.table("~/selection/analysis/poly/Height/pred_height.txt", as.is=TRUE)
data[,1] <- gsub("_", " ", data[,1])
pops <- unique(data[,1])

cols <- brewer.pal( length(pops), "Set1")
cols[6] <- "darkgrey"

pdf(paste0(out, outname), width=12, height=6)
plot(0,0, col="white", xlim=c(-8000, 0), ylim=ylim, bty="n", xlab="Years before present", ylab="Genetic height", xaxt="n")

for(i in 1:length(pops)){
    pop <- pops[i]
    
    int.starts <- data[data[,1]==pop,3]
    int.ends <- data[data[,1]==pop,4]
    ht <- data[data[,1]==pop,5]
    rect(-int.starts, ht-0.02, -int.ends, ht+0.02, col=cols[i], density=20, angle=ang*i, border=cols[i])
    ln <- length(int.ends)

    if(ln>1){
        segments(-int.ends[1:(ln-1)], ht[1:(ln-1)], -int.starts[2:ln], ht[2:ln], col=cols[i], lwd=2)
    }
}
legend("bottomright", pops, col=cols,lwd=2, bty="n")
axis(1, at=-seq(8,0,-2)*1000, labels=format(seq(8,0,-2)*1000, big.mark=","))

dev.off()
