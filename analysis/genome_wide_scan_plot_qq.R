#Make a Manhattan plot of the scan results
source("~/selection/code/lib/mh_plot_lib.R")
library(RColorBrewer)

##############################################################

results <- "~/selection/analysis/gscan/scan_results.txt"

##############################################################

results <- read.table(results, header=TRUE)

data <- read.table("~/data/v6/use/v61kg_europe2names.snp", as.is=TRUE)

rownames(results) <- results$ID

res <- data.frame(ID=results$ID, PVAL=results$corrected.p)
dat <- data[,c(1,2,4)]
colnames(dat) <- c("ID", "CHR", "POS")
res <- merge(res,dat,by="ID")
res <- res[order(res$CHR, res$POS),]

selection <- read.table("~/selection/data/cms/All_selection.snp", as.is=TRUE)
cms.EUR <- read.table("~/selection/data/cms/CMS_Tabrizi_CEU.snp", as.is=TRUE)
res$TAG <- paste(res$CHR, res$POS, sep=":")
selection$TAG <- paste(selection[,2], selection[,4], sep=":")
cms.EUR$TAG <- paste(cms.EUR[,2], cms.EUR[,4], sep=":")

#First show the effect of genomic control on "neutral" snps.
neutral <- res[!(res$TAG %in% selection$TAG),]
selection <- res[(res$TAG %in% selection$TAG)&!(res$TAG %in% cms.EUR$TAG),]

png(paste0("~/selection/analysis/gscan/qq_plot_all.png"), width=400, height=400)
par(mar=c(4,4,1,1))
qqPlotOfPValues(neutral$PVAL, col="#377EBA", linecol="black", ci=FALSE)
draw.qq.ci(NROW(neutral), border.col="#377EBA20", fill.col="#377EBA20")
qqPlotOfPValues(selection$PVAL, col="#E41A1C", add=T)
draw.qq.ci(NROW(selection), border.col="#E41A1C20", fill.col="#E41A1C20")
dev.off()

## png(paste0("~/selection/gscan/qq_zoom", tag ,".png"), width=600, height=600)
## qqPlotOfPValues(neutral$PVAL, col="#377EBA", linecol="black", xlim=c(0,5), ylim=c(0,5), ci=FALSE)
## draw.qq.ci(NROW(neutral), border.col="#377EBA80", fill.col="#377EBA20")
## qqPlotOfPValues(selection$PVAL, col="#E41A1C", add=T)
## draw.qq.ci(NROW(selection), border.col="#E41A1C80", fill.col="#E41A1C20")
## dev.off()

#Now plot the selection snps spit up by functional category.
cat.data <- read.table("~/selection/data/Selection_snps_by_category.txt", as.is=TRUE, header=TRUE)
cat.data$GWAS <- ifelse(cat.data$GWAS|cat.data$Pheno1|cat.data$Pheno2,1,0)
cat.data$Immune <- cat.data$Pheno3
cat.data$TAG <- paste(cat.data$CHR, cat.data$POS, sep=":")

png(paste0("~/selection/analysis/gscan/qq_plot_cat.png"), width=400, height=400)
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
