#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
source("~/selection/code/lib/3pop_lib.R")

#Modern GBR, CEU, IBS, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
root <- "~/selection/counts/all"
out <- "~/selection/analysis/gscan/"
results.tag <- ""

lambda=NA

########################################################################

#Compute the likelihood 

pops <- c("WHG", "LBK_EN", "Yamnaya", "CEU", "GBR", "IBS", "TSI")
## A <- matrix(c(0.187, 0.312, 0.501, 0.160, 0.413, 0.427, 0, 0.764, 0.236, 0, 0.714, 0.286),3, 4)
A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.226, 0, 0.712, 0.287),3, 4) 

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)

## Merge WHG (not inbred so counted separately...)
counts$WHG <- counts$Loschbour+counts$LaBrana1+counts$HungaryGamba_HG
totals$WHG <- totals$Loschbour+totals$LaBrana1+totals$HungaryGamba_HG

## counts$LBK_EN <- counts$LBK_EN +counts$Stuttgart+counts$HungaryGamba_EN+counts$Spain_EN
## totals$LBK_EN <- totals$LBK_EN +totals$Stuttgart+totals$HungaryGamba_EN+totals$Spain_EN

## cat("Using Falchenstein!!\n")
## results.tag <- "_with_Falchenstein"
## counts$WHG <- counts$Loschbour+counts$LatvianMesolithic+counts$LaBrana1+counts$HungaryGamba_HG+counts$GermanMesolithic
## totals$WHG <- totals$Loschbour+totals$LatvianMesolithic+totals$LaBrana1+totals$HungaryGamba_HG+totals$GermanMesolithic

cat("Using all the EN!!\n")
counts$LBK_EN <- counts$LBK_EN +counts$Stuttgart+counts$HungaryGamba_EN+counts$Spain_EN+counts$Starcevo_EN+counts$LBKT_EN
totals$LBK_EN <- totals$LBK_EN +totals$Stuttgart+totals$HungaryGamba_EN+totals$Spain_EN+totals$Starcevo_EN+totals$LBKT_EN

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
    
    results[i,] <- test.3pop(N, N.A, A)
}

results <- results[!is.na(results[,2]),]
colnames(results) <- c("ChiSq", "uncorrected.p")
results <- data.frame(results)

#Correct according to neutral snps only
selection <- read.table("~/selection/data/cms/All_selection.snp", as.is=TRUE)
neutral <- !(rownames(results) %in% selection[,1])

png(paste0("~/selection/analysis/gscan/scan_gc", results.tag ,".png"), width=800, height=800)
degf <- dim(A)[2]
qq.exp.pts <- rexp(NROW(results))/log(10)
qqplot(qq.exp.pts, -log10(results[,2]), pch=1, col="red", xlab="Expected", ylab="Observed", bty="n", xlim=c(0,ceiling(log10(NROW(data)))), cex=0.5)
abline(0,1,col="black")
if(is.na(lambda)){
    lambda <- median(results[neutral,"ChiSq"])/qchisq(0.5, df=degf)
    cat(paste0("lambda=", lambda))
}
corrected.p <- pchisq(results[,1]/lambda, df=degf, lower.tail=F)
tt <- qqplot(qq.exp.pts, -log10(corrected.p), plot.it=FALSE)
points(tt, pch=16, col="darkblue", bty="n")
dev.off()

results <- cbind(rownames(results), results, corrected.p)
colnames(results) <- c("ID", "ChiSq", "uncorrected.p", "corrected.p")
results <- data.frame(results)
write.table(results, paste0("~/selection/analysis/gscan/scan_results", results.tag, ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

