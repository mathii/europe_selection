#First run genome_wide_scan.R to generate the results
source("~/selection/code/lib/mh_plot_lib.R")

##############################################################

results.tag <- ""
results <- paste0("~/selection/analysis/gscan/scan_results", results.tag, ".txt")

##############################################################

results <- read.table(results, header=TRUE)

## QQ plots

## Manhattan plot
data <- read.table("~/data/v6/use/v61kg_europe2names.snp", as.is=TRUE)
dist <- 100000

res <- data.frame(ID=results$ID, PVAL=results$corrected.p)
dat <- data[,c(1,2,4)]
colnames(dat) <- c("ID", "CHR", "POS")
res <- merge(res,dat,by="ID")
res <- res[order(res$CHR, res$POS),]


## Manhattan plot
png(paste0("~/selection/analysis/gscan/mh_plot", results.tag ,".png"), width=800, height=400)
par(mar=c(2,4,1,1))
MH.plot(res, color.loci=data.frame())
abline(h=6.79, col="red", lty=2)
dev.off()
