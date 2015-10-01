## box series plot of height. 
## Just going to hardcode the results in here for now.

library(plotrix)
library(RColorBrewer)

version <- NA
cA <- commandArgs(TRUE)
if(length(cA)){
  version <- cA[1]
  snplist <- cA[2]
}

########################################################################
## Details
out <-  paste0("~/selection/analysis/", version, "/series/heightseries.",snplist,".pdf")
height.values <-  paste0("~/selection/analysis/", version, "/series/height_series_mcmc_estimates.",snplist,".txt")
dates <- paste0("~/selection/code/files/",version,"/population_dates.txt")
colfile <- paste0("~/selection/code/files/",version,"/population_cols.txt")

ang <- 20

########################################################################
## Details
dates <- read.table(dates, as.is=TRUE)
colnames(dates) <- c("MPOP", "POP", "START", "END")
heights <- read.table(height.values, as.is=TRUE, header=TRUE)
heights$POP <- rownames(heights)
data <- merge(dates, heights, by="POP")

data[,1] <- gsub("_", " ", data[,"MPOP"])
pops <- unique(data[,"MPOP"])
data <- data[order(data$MPOP, data$END, data$START),]

if(file.exists(colfile)){
  cc <- read.table(colfile, as.is=TRUE, comment.char="")
  cols <- cc[,2]
  names(cols) <- cc[,1]
} else{
  if(length(pops)<=9){
    cols <- brewer.pal( length(pops), "Set1")
    cols[6] <- "darkgrey"
  } else{
    cols <- c(brewer.pal( length(pops), "Set1"), brewer.pal( length(pops), "Set2"))
    cols[6] <- "darkgrey"
  }
  names(cols) <- pops
}

a.i <- c(data[,"Post.5"], data[,"Post.95"], data[,"MLE"])
ylim <- c(floor(min(a.i)*10)/10, max(a.i)*10/10)


pdf(out, width=12, height=6)
plot(0,0, col="white", xlim=c(-8500, 0), ylim=ylim, bty="n", xlab="Years before present", ylab="Genetic height", xaxt="n")

for(i in 1:length(pops)){
    pop <- pops[i]

    kk=data[,"MPOP"]==pop
    int.starts <- data[kk,"END"]
    int.ends <- data[kk,"START"]
    ht <- data[kk,"Post.Mn"]
    ## rect(-int.starts, ht-0.02, -int.ends, ht+0.02, col=cols[i], density=20, angle=ang*i, border=cols[i])
    rect(-int.starts, data[kk,"Post.5"], -int.ends, data[kk,"Post.95"], col=cols[pop], density=20, angle=ang*i, border=cols[pop])
    text(-int.ends, data[kk,"Post.95"], data[kk,"POP"], adj=c(0,-0.2))
    points(-0.5*(int.starts+int.ends), data[kk,"MLE"], cex=2, pch=16, col=cols[pop])
                
    ## ln <- length(int.ends)
    ## if(ln>1 & pop!="Hunter_Gatherers"){
    ##     segments(-int.ends[1:(ln-1)], ht[1:(ln-1)], -int.starts[2:ln], ht[2:ln], col=cols[i], lwd=2)
  ## }
}
## legend("bottomright", pops, col=cols[pops],lwd=2, bty="n")
axis(1, at=-seq(8,0,-2)*1000, labels=format(seq(8,0,-2)*1000, big.mark=","))

dev.off()
