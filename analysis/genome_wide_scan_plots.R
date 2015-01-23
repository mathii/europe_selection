#Manhattan and qq plots
#First run genome_wide_scan{*}.R to generate the results
library(RColorBrewer)

source("~/selection/code/lib/mh_plot_lib.R")

##############################################################

results.tag <- ""
if(length(commandArgs(TRUE))){
    results.tag <- commandArgs(TRUE)[1]
}

results <- paste0("~/selection/analysis/gscan/scan_results", results.tag, ".txt")
snpdata <- "~/data/v6/use/v61kg_europe2names.snp"
degf <- 4

##############################################################

## Read inputs
results <- read.table(results, header=TRUE, as.is=TRUE)
selection <- read.table("~/selection/data/cms/All_selection.snp", as.is=TRUE)
data <- read.table(snpdata, as.is=TRUE)

## Apply genomic control if we haven't already
if(!("corrected.p"  %in% names(results))){
    neutral <- !(results[,"ID"] %in% selection[,1])
    lambda <- median(results[neutral,"ChiSq"])/qchisq(0.5, df=degf)
    corrected.p <- pchisq(results[,"ChiSq"]/lambda, df=degf, lower.tail=F)
    results <- cbind(results, corrected.p)
    write.table(results, paste0("~/selection/analysis/gscan/scan_results", results.tag, ".txt"),
                row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    cat(paste0("Used lambda = ", lambda, "\n"))
}else{cat("Using corrected p-values\n")}


## Merge SNP position data with results
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

## QQ plots
res$TAG <- paste(res$CHR, res$POS, sep=":")
selection$TAG <- paste(selection[,2], selection[,4], sep=":")

#First show the effect of genomic control on "neutral" snps.
neutral <- res[!(res$TAG %in% selection$TAG),]
selection <- res[(res$TAG %in% selection$TAG),]
png(paste0("~/selection/analysis/gscan/qq_plot_all", results.tag, ".png"), width=400, height=400)
par(mar=c(4,4,1,1))
qqPlotOfPValues(neutral$PVAL, col="#377EBA", linecol="black", ci=FALSE)
draw.qq.ci(NROW(neutral), border.col="#377EBA20", fill.col="#377EBA20")
qqPlotOfPValues(selection$PVAL, col="#E41A1C", add=T)
draw.qq.ci(NROW(selection), border.col="#E41A1C20", fill.col="#E41A1C20")
dev.off()

#Now plot the selection snps spit up by functional category.
cat.data <- read.table("~/selection/data/Selection_snps_by_category.txt", as.is=TRUE, header=TRUE)
cat.data$GWAS <- ifelse(cat.data$GWAS|cat.data$Pheno1|cat.data$Pheno2,1,0)
cat.data$Immune <- cat.data$Pheno3
cat.data$TAG <- paste(cat.data$CHR, cat.data$POS, sep=":")

png(paste0("~/selection/analysis/gscan/qq_plot_cat", results.tag, ".png"), width=400, height=400)
par(mar=c(4.1,4.1,1,1))
cats <- c("GWAS", "CMS", "HiDiff", "Immune", "HLA")
cols <- brewer.pal(length(cats), "Set1")
for(i in 1:length(cats)){
    this.lot <- res[res$TAG %in% cat.data$TAG[cat.data[,cats[i]]==1],]
    qqPlotOfPValues(this.lot$PVAL, col=cols[i], ci=FALSE, add=(i!=1), xlim=c(0,5), ylim=c(0,5),  linecol="black" )
    draw.qq.ci(NROW(this.lot), border.col=paste0(cols[i],"20"), fill.col=paste0(cols[i],"20"))
}
legend("topleft", cats, col=cols, pch=16, bty="n")
dev.off()
