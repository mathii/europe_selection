#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
source("~/selection/code/lib/3pop_lib.R")
source("~/selection/code/lib/mh_plot_lib.R")

#Modern GBR, CEU, FIN, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
root <- "~/selection/counts/all"
out <- "~/selection/analysis/gscan/"

lambda=NA

########################################################################

#Compute the likelihood 

pops <- c( "CEU", "GBR", "IBS", "TSI")

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)

data <- counts[,1:5]
counts <- data.matrix(counts[,6:NCOL(counts)])
totals <- data.matrix(totals[,6:NCOL(totals)])
counts <- counts[,pops]
totals <- totals[,pops]

results <- matrix(0, nrow=NROW(data), ncol=2)
rownames(results) <- data$ID

for(i in 1:NROW(data)){
    cat(paste0("\r", i))
    N <- totals[i,]
    N.A <- counts[i,]

    if(all(N==N.A)){
        results[i,] <- NA
        next
    }
    
    results[i,] <- test.diff(N, N.A)
}

results <- results[!is.na(results[,2]),]

degf <- length(N)-1
qqplot(rexp(NROW(results))/log(10), -log10(results[,2]), pch=1, col="red", xlab="Expected", ylab="Observed", bty="n", xlim=c(0,ceiling(log10(NROW(data)))), cex=0.5)
abline(0,1,col="black")
if(is.na(lambda)){
    lambda <- median(results[,1])/qchisq(0.5, df=degf)
    cat(paste0("lambda=", lambda))
}
corrected.p <- pchisq(results[,1]/lambda, df=degf, lower.tail=F)
tt <- qqplot(rexp(NROW(results))/log(10), -log10(corrected.p), plot.it=FALSE)
points(tt, pch=16, col="darkblue", bty="n")

results <- cbind(rownames(results), results, corrected.p)
colnames(results) <- c("ID", "ChiSq", "uncorrected.p", "corrected.p")
results <- data.frame(results)
write.table(results, "~/selection/analysis/gscan/scan_diff_results.txt", row.names=T, col.names=T, quote=F)
