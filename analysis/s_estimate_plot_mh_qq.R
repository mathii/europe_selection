#Make a manhattan plot from the s_estimates.
source("~/selection/code/lib/mh_plot_lib.R")

results <- read.table("~/selection/analysis/s_estimates/ALL_s_estimates.txt", as.is=TRUE, header=FALSE)
colnames(results) <- c("ID", "s.est", "p", "ci.lower", "ci.upper")
data <- read.table("~/data/v6/use/v61kg_europe2names.snp", as.is=TRUE)
data <- data[,c(1,2,4)]
colnames(data) <- c("ID", "CHR", "POS")
results.tag <- ""

png(paste0("~/selection/analysis/s_estimates/qq_plot", results.tag ,".png"), width=800, height=400)
qq.exp.pts <- rexp(NROW(results))/log(10)
qqplot(qq.exp.pts, -log10(results$p), pch=1, col="red", xlab="Expected", ylab="Observed", bty="n", cex=0.5)
abline(0,1,col="black")
dev.off()

res <- data.frame(ID=results$ID, PVAL=results$p)
res <- merge(res,data,by="ID")
res <- res[order(res$CHR, res$POS),]
res <- res[!is.na(res$PVAL),]

png(paste0("~/selection/analysis/s_estimates/mh_plot", results.tag ,".png"), width=800, height=400)
par(mar=c(2,4,1,1))
MH.plot(res)
abline(h=6.79, col="red", lty=2)

dev.off()
