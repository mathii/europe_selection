#Manhattan and qq plots
#First run genome_wide_scan{*}.R to generate the results
library(RColorBrewer)

source("~/selection/code/lib/mh_plot_lib.R")

##############################################################

results.tag <- ""
version <- "" 
degf <- NA                      #degrees of freedom for test
what <- "gscan"                #plot other things than gscan
cA <- commandArgs(TRUE)
cutoff <- 0
if(length(cA)){
    results.tag <- cA[1]
    version <- cA[2]
    degf <- as.numeric(cA[3])
    if(length(cA)>3){what <- cA[4]}
}else{
    stop("Must specify results tag as first argument")
}

if(version==""){stop("Must specify version as second argument")}
if(is.na(degf)){stop("Must specify degrees of freedom as third argument")}

results <- paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag, ".txt.gz")
snpdata <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.snp")

logfile <- paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag, ".log")
cat("SCRIPT: genome_wide_scan_plots.R\n", file=logfile)
cat(paste0(paste("ARGS:", paste(cA)),"\n"), file=logfile, append=TRUE)

##############################################################

## Read inputs
results <- read.table(results, header=TRUE, as.is=TRUE)

selection <- read.table(paste0("~/data/",version, "/use/All_selection.snp"), as.is=TRUE)
neutral <- !(results[,"ID"] %in% selection[,1])
data <- read.table(snpdata, as.is=TRUE)

include <- results$ChiSq>0
cat(paste0("IGNORE: ", sum(!include), " SNPs with non-positive test statistics\n"), file=logfile, append=TRUE)

## Apply genomic control if we haven't already
## if(!("corrected.p"  %in% names(results))){
lambda <- median(results[neutral&include,"ChiSq"])/qchisq(0.5, df=degf)
corrected.p <- pchisq(results[,"ChiSq"]/lambda, df=degf, lower.tail=F)
results <- cbind(results, corrected.p)
## write.table(results, paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag, ".txt"),
##             row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
cat(paste0("LAMBDA: ", lambda, "\n"), file=logfile, append=TRUE)
## }else{cat("Using corrected p-values\n")}

## Merge SNP position data with results
res <- data.frame(ID=results$ID, PVAL=results$corrected.p, stringsAsFactors=FALSE)
res <- res[include,]
dat <- data[,c(1,2,4)]
colnames(dat) <- c("ID", "CHR", "POS")
res <- merge(res,dat,by="ID")
res <- res[order(res$CHR, res$POS),]

sig.level <- -log10(0.05/NROW(res))
lo.sig <- sig.level-2
cat(paste0("GWSIG: ", sig.level, "\n"),file=logfile, append=TRUE)
cat(paste0("LOSIG: ", lo.sig, "\n"),file=logfile, append=TRUE)
if(version=="v8"){
    mixmap <- read.table(paste0("~/selection/code/files/v8/mixtures/Choice", gsub( "_read", "", results.tag)), as.is=TRUE, header=FALSE)
    pops <- unique(mixmap[,2])
    for(i in 1:length(pops)){
        cat(paste0("MPOP", i, ": ", paste(mixmap[mixmap[,2]==pops[i],1], collapse=","), "\n"),file=logfile, append=TRUE)
    }
}                       
if(all(c("eff.N1", "eff.N2", "eff.N3") %in% names(results))){
    cat(paste0("EFSIZ: ", paste(round(colMeans(results[,c("eff.N1", "eff.N2", "eff.N3")]),2), collapse=" "), "\n"),file=logfile, append=TRUE)
}

## Manhattan plot
png(paste0("~/selection/analysis/",version,"/", what ,"/mh_plot", results.tag ,".png"), width=800, height=400)
par(mar=c(2,4,1,1))
MH.plot(res, color.loci=data.frame())
abline(h=sig.level, col="red", lty=2)
dev.off()

## Per-chromosome Manhattan plots
sig.chrs <- unique(res$CHR[res$PVAL<(10^-(sig.level))])
for(c in sig.chrs){
    png(paste0("~/selection/analysis/",version,"/", what ,"/mh_plot", results.tag ,".chr", c, ".png"), width=250, height=250)
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
cat.data <- read.table(paste0("~/data/",version, "/use/Cat_selection.snp"), as.is=TRUE, header=TRUE)
cat.data$GWAS <- ifelse(cat.data$GWAS|cat.data$Pheno1|cat.data$Pheno2,1,0)
cat.data$Immune <- cat.data$Pheno3
cat.data$TAG <- paste(cat.data$CHR, cat.data$POS, sep=":")

png(paste0("~/selection/analysis/",version,"/", what ,"/qq_plot_cat", results.tag, ".png"), width=400, height=400)
par(mar=c(4.1,4.1,1,1))
cats <- c("GWAS", "CMS", "HiDiff", "Immune", "HLA", "eQTL")
cats <- cats[cats %in% colnames(cat.data)]
cols <- brewer.pal(length(cats), "Set1")
for(i in 1:length(cats)){
    this.lot <- res[res$TAG %in% cat.data$TAG[cat.data[,cats[i]]==1],]
    qqPlotOfPValues(this.lot$PVAL, col=cols[i], ci=FALSE, add=(i!=1), xlim=c(0,4), ylim=c(0,20),  linecol="black" )
    draw.qq.ci(NROW(this.lot), border.col=paste0(cols[i],"20"), fill.col=paste0(cols[i],"20"))
}
legend("topleft", cats, col=cols, pch=16, bty="n")
dev.off()

isig <- indep.signals(res, 10^(-lo.sig), 10^(-sig.level), 1e6)

if(file.exists("~/selection/data/genes/refseq_inm.txt")){
    gene.names<-read.table("~/selection/data/genes/refseq_inm.txt", header=TRUE, as.is=TRUE)
    isig <- annotate.with.genes(isig, gene.names)
}
   
write.table(isig, paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag ,".signals.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(isig[isig$n.sig>2|isig$n.gw.sig>1,], paste0("~/selection/analysis/",version,"/", what ,"/scan_results", results.tag ,".clean_signals.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## Now make a cleaned Manhatten plot, remiving everything that's genome-wide significant
## But not supported by anything within two p-value orders of magnitude. 
to.remove <- isig$lead.snp[isig$n.sig<2]
cat(paste0("REMOVE: ", sum(isig$n.sig<2), "\n"),file=logfile, append=TRUE)
clean.res <- res
clean.res <- clean.res[!(clean.res$ID %in% to.remove),]
png(paste0("~/selection/analysis/",version,"/", what ,"/mh_plot", results.tag ,".cleaned.png"), width=800, height=400)
par(mar=c(2,4,1,1))
MH.plot(clean.res, color.loci=data.frame())
abline(h=sig.level, col="red", lty=2)
csig <- isig[isig$n.sig>2 & isig$n.gw.sig>1,]

dev.off()




