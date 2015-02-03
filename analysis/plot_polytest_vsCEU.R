#Plot the results of the polygenic selection test against CEU.
## library(shape)
library(RColorBrewer)

result.tag <- ""
if(length(commandArgs(TRUE))){
    result.tag <- commandArgs(TRUE)[1]
}

traits <- c("Height", "BMI", "WHR", "T2D", "IBD", "Lipids")

polypops <- scan("~/selection/data/polypops.txt", what="")
popkey <- letters[1:length(polypops)]
names(popkey) <- polypops

xlim=c(-6,4)
xticks <- c(-6, -4, -2, 0, 2, 4)

pdf(paste0("~/selection/analysis/poly/results/vsCEU", result.tag ,".pdf"), height=5, width=6)
par(mar=c(4.1,4.1,4.1,2.1))
plot(0,0, col="white", bty="n", xaxt="n", yaxt="n", bty="n", xlim=xlim+c(0,2), ylim=c(-0.5, length(traits)+0.5), xlab="", ylab="")

i=length(traits)
for(trait in traits){
    results <- read.table(paste0("~/selection/analysis/poly/", trait, "/", trait, "_vsCEU", result.tag , ".results.txt"), as.is=T, header=T)
    pops <- sapply(strsplit(results[,1],","), "[[", 2)
    Z <- sqrt(results[,2])
    gv <- matrix(as.numeric(unlist(strsplit(results[,5], ","))), ncol=2, byrow=TRUE)
    Z <- Z*ifelse(gv[,2]>gv[,1],1,-1)
#    Arrows(xlim[1], i, xlim[2], i, code=3)
    segments(xlim[1], i, xlim[2], i, lwd=)
    thecols <- ifelse(results[,3]<0.01, "Red", ifelse(results[,3]<0.05, "Purple", "Blue"))
    points(Z, i+0*Z, pch=21, bg=thecols, col="grey", cex=1.5)

    lab=ifelse(results[,3]<0.01, pops, "")

    zz <- order(Z)
    offset <- -1+2*cumsum(abs(Z[zz])>3)%%2
    if(offset[1]<0){offset <- 2-offset}
    text(Z[zz], i-0.02, lab[zz], pos=ifelse(offset>0,1,3), cex=0.8)
    text(xlim[2]-0.2, i, trait, pos=4, cex=1.2)
    i=i-1
}

segments(xlim[1], 0, xlim[2],0)
xlabs <- ifelse(xticks<=0, as.character(xticks), paste0( "+", xticks))
for(xx in xticks){
    abline(v=xx, col="#00000020")
}
text(xticks, -0.1, xlabs, pos=1)
dev.off()
