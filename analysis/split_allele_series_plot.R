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
pp[,1]<-gsub("_", " ", pp[,1])
########################################################################
## Details
snpname <- "rs12913832"
flip <- TRUE
outname <- "herc2series.pdf"

data <- read.table("~/selection/counts/all.reads.freq", as.is=TRUE, header=TRUE)
freq <- data[,6:NCOL(data)]
data <- data[,1:5]
rownames(freq) <- data[,1]

uci <- read.table("~/selection/counts/all.reads.highCI.freq", as.is=TRUE, header=TRUE)
uci <- unlist(uci[uci[,1]==snpname,6:NCOL(uci)])
lci <- read.table("~/selection/counts/all.reads.lowCI.freq", as.is=TRUE, header=TRUE)
lci <- unlist(lci[lci[,1]==snpname,6:NCOL(lci)])


pops <- unique(pp[,1])

ff<-freq[snpname,]
if(flip){
    ff <- 1-ff
    tmp1 <- 1-uci
    uci <- 1-lci
    lci <- tmp1
}
pp$f <- unlist(ff[pp[,2]])
pp$lci <- unlist(lci[pp[,2]])
pp$uci <- unlist(uci[pp[,2]])
         
cols <- brewer.pal( length(pops), "Set1")
cols[6] <- "darkgrey"

pdf(paste0(out, outname), width=12, height=6)
plot(0,0, col="white", xlim=c(-8000, 0), ylim=ylim, bty="n", xlab="Years before present", ylab="HERC2 allele frequency", xaxt="n")

for(i in 1:length(pops)){
    pop <- pops[i]
    
    int.starts <- pp[pp[,1]==pop,3]
    int.ends <- pp[pp[,1]==pop,4]
    int.mids <- 0.5*(int.starts+int.ends)
    if(i==3){int.mids[1] <- int.mids[1]+20}
    af <- pp[pp[,1]==pop,"f"]
    uc <- pp[pp[,1]==pop,"uci"]
    lc <- pp[pp[,1]==pop,"lci"]

    rect(-int.starts, af-0.02, -int.ends, af+0.02, col=cols[i], density=20, angle=ang*i, border=cols[i])
    ln <- length(int.ends)

    inc.up <- af<0.98
    if(sum(inc.up)){
        segments(-int.mids[inc.up], af[inc.up]+0.02, -int.mids[inc.up], uc[inc.up], col=cols[i], lwd=2)
    }
    inc.dn <- af>0.02
    if(sum(inc.dn)){
        segments(-int.mids[inc.dn], af[inc.dn]-0.02, -int.mids[inc.dn], lc[inc.dn], col=cols[i], lwd=2)
    }
    
    if(ln>1){
        segments(-int.ends[1:(ln-1)], af[1:(ln-1)], -int.starts[2:ln], af[2:ln], col=cols[i], lwd=2)
    }
}
legend("bottomright", pops, col=cols,lwd=2, bty="n")
axis(1, at=-seq(8,0,-2)*1000, labels=format(seq(8,0,-2)*1000, big.mark=","))
dev.off()
