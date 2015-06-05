#Manhattan and qq plots
#First run genome_wide_scan{*}.R to generate the results
library(RColorBrewer)

source("~/selection/code/lib/mh_plot_lib.R")

##############################################################

results.tag <- ""
version <- "" 
degf <- 4
what <- "gsan"
cA <- commandArgs(TRUE)
if(length(cA)){
    results.tag <- cA[1]
    version <- cA[2]
    degf <- as.numeric(cA[3])
    if(length(cA)>3){what <- cA[4]}
}

results <- paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag, ".txt")
snpdata <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.snp")

##############################################################

## Read inputs
results <- read.table(results, header=TRUE, as.is=TRUE)
selection <- read.table("~/selection/data/cms/All_selection.snp", as.is=TRUE)
neutral <- !(results[,"ID"] %in% selection[,1])
data <- read.table(snpdata, as.is=TRUE)

include <- results$ChiSq>0
cat(paste0("Ignoring ", sum(!include), " SNPs with non-positiveq test statistics"))

## Apply genomic control if we haven't already
if(!("corrected.p"  %in% names(results))){
    lambda <- median(results[neutral&include,"ChiSq"])/qchisq(0.5, df=degf)
    corrected.p <- pchisq(results[,"ChiSq"]/lambda, df=degf, lower.tail=F)
    results <- cbind(results, corrected.p)
    write.table(results, paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag, ".txt"),
                row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    cat(paste0("Used lambda = ", lambda, "\n"))
}else{cat("Using corrected p-values\n")}

## Merge SNP position data with results
res <- data.frame(ID=results$ID, PVAL=results$corrected.p)
res <- res[include,]
dat <- data[,c(1,2,4)]
colnames(dat) <- c("ID", "CHR", "POS")
res <- merge(res,dat,by="ID")
res <- res[order(res$CHR, res$POS),]


sig.level <- -log10(0.05/NROW(res))
cat(paste0("Used sig.level = ", sig.level, "\n"))

## Manhattan plot
## png(paste0("~/Dropbox/posters/AAPA/mh_plot", results.tag ,".png"), width=16, height=8, units="in", res=400)
## par(bg = "grey92")
png(paste0("~/selection/analysis/",version,"/", what ,"/mh_plot", results.tag ,".png"), width=800, height=400)
par(mar=c(2,4,1,1))
MH.plot(res, color.loci=data.frame())
abline(h=sig.level, col="red", lty=2)
dev.off()

## Per-chromosome Manhattan plots
sig.chrs <- unique(res$CHR[res$PVAL<(10^-(sig.level))])
for(c in sig.chrs){
    png(paste0("~/selection/analysis/",version,"/", what ,"/mh_plot", results.tag ,".chr", c, ".png"), width=250, height=250)
    ## png(paste0("~/Dropbox/posters/AAPA/mh_plot.chr", c ,".png"), width=4, height=4, units="in", res=400)
    ## par(bg = "grey92")
    par(mar=c(1,2,1,1))
    MH.plot(res[res$CHR==c,], color.loci=data.frame(), chr.labels=FALSE)
    abline(h=sig.level, col="red", lty=2)
    dev.off()
}

## QQ plots
res$TAG <- paste(res$CHR, res$POS, sep=":")
selection$TAG <- paste(selection[,2], selection[,4], sep=":")

#First show the effect of genomic control on "neutral" snps.
neutral <- res[!(res$TAG %in% selection$TAG),]
selection <- res[(res$TAG %in% selection$TAG),]
png(paste0("~/selection/analysis/",version,"/", what ,"/qq_plot_all", results.tag, ".png"), width=400, height=400)
par(mar=c(4,4,1,1))
qqPlotOfPValues(neutral$PVAL, col="#377EBA", linecol="black", ci=FALSE, xlim=c(0,6), ylim=c(0,30))
draw.qq.ci(NROW(neutral), border.col="#377EBA20", fill.col="#377EBA20")
qqPlotOfPValues(selection$PVAL, col="#E41A1C", add=T)
draw.qq.ci(NROW(selection), border.col="#E41A1C20", fill.col="#E41A1C20")
dev.off()

#Now plot the selection snps spit up by functional category.
cat.data <- read.table("~/selection/data/Selection_snps_by_category.txt", as.is=TRUE, header=TRUE)
cat.data$GWAS <- ifelse(cat.data$GWAS|cat.data$Pheno1|cat.data$Pheno2,1,0)
cat.data$Immune <- cat.data$Pheno3
cat.data$TAG <- paste(cat.data$CHR, cat.data$POS, sep=":")

png(paste0("~/selection/analysis/",version,"/", what ,"/qq_plot_cat", results.tag, ".png"), width=400, height=400)
par(mar=c(4.1,4.1,1,1))
cats <- c("GWAS", "CMS", "HiDiff", "Immune", "HLA")
cols <- brewer.pal(length(cats), "Set1")
for(i in 1:length(cats)){
    this.lot <- res[res$TAG %in% cat.data$TAG[cat.data[,cats[i]]==1],]
    qqPlotOfPValues(this.lot$PVAL, col=cols[i], ci=FALSE, add=(i!=1), xlim=c(0,4), ylim=c(0,20),  linecol="black" )
    draw.qq.ci(NROW(this.lot), border.col=paste0(cols[i],"20"), fill.col=paste0(cols[i],"20"))
}
legend("topleft", cats, col=cols, pch=16, bty="n")
dev.off()
