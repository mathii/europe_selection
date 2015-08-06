## Merge the plots of robustness to admixture and Proportion

version <- "v8"
what <- "Admixture"
seeds <- seq(101,125)

results <- list()
for(i in 1:length(seeds)){
    results[[i]] <- read.table(paste0("~/selection/analysis/v8/power/reads_robust_power_", what,"_seed_", seeds[i],"_lambda.txt"), as.is=TRUE, header=TRUE)
}

result <- results[[1]]
for(i in 2:length(seeds)){
    result <- result +results[[i]]
}
result <- result/length(seeds)

write.table(result, paste0("~/selection/analysis/v8/power/reads_robust_power_", what, ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

pdf(paste0("~/selection/analysis/", version,"/power/reads_robust_power_",what,".pdf"))
par(mar=c(5,4,4,4))
plot(result[,1], result[,2], col="#377EBA", type="b", pch=16, bty="n", lwd=3, xlab="Random proportion", ylab="Genomic inflation factor", yaxt="n", xaxt="n", ylim=c(1.35,1.6))
axis(1, lwd=2)
axis(2, col="#377EBA", lwd=2)
par(new=TRUE)
plot(result[,1], result[,3], col="#CC5500", type="b", pch=16, bty="n", lwd=3, axes=FALSE, xlab="", ylab="", ylim=c(0.15, 0.3))
axis(4, col="#CC5500", lwd=2)
mtext("Power", 4, line=3)
dev.off()
