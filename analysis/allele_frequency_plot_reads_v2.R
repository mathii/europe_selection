## allele frequency plot that plots all the ancient population frequencies, 
## along with the range of the modern population frequecies.

source("~/selection/code/lib/readlib.R")
source("~/selection/code/lib/plotlib.R")
source("~/selection/code/lib/3pop_lib.R")

library(plotrix)
library(RColorBrewer)

########################################################################
## Details
ylim <- c(0,1)
outsize=c(6,6)
error.prob <- 0.001

version="v6"
what <- "figure2"
map <- "2"

cA <- commandArgs(TRUE)
if(length(cA)>0){version <- cA[1]}
if(length(cA)>1){what <- cA[2]}
if(length(cA)>1){map <- cA[3]}

outname <- paste0(what, "map", map, "v2.pdf")
readmefile <- paste0("~/selection/code/files/",version, "/",  what, ".readme")
root <- paste0("~/selection/counts/",version,"/all")
read.root <- paste0("~/data/",version,"/reads/")
indfile <- paste0("~/data/",version,"/use/",version,"1kg_europe2names.ind")
out <- paste0("~/selection/analysis/",version,"/series/")

## ########################################################################

if(version=="v6"){
int.names <- c("SHG", "WHG", "EN", "MN", "SteppeBA","LN/BA")
long.names <- c("SHG", "WHG", "ENeo", "MNeo", "SteppeBA", "LNeo/BA")
int.include <- c("SHG", "WHG", "WHG", "WHG", "EN", "EN", "EN", "EN", "EN", "EN", "MN", "MN", "MN", "MN", "SteppeBA", "SteppeBA", "SteppeBA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "LN/BA", "CEU/IBS", "CEU/IBS")
names(int.include) <-c("Motala_HG", "Loschbour", "Iberian_Mesolithic", "HungaryGamba_HG", "Starcevo_EN", "Stuttgart", "Spain_EN", "LBK_EN", "LBKT_EN", "HungaryGamba_EN", "Spain_MN", "Baalberge_MN", "Iceman", "Esperstedt_MN", "Yamnaya", "Poltavka", "Srubnaya", "HungaryGamba_CA", "Alberstedt_LN", "Corded_Ware_LN", "Bell_Beaker_LN", "BenzigerodeHeimburg_LN", "Unetice_EBA", "HungaryGamba_BA", "Halberstadt_LBA")
include.reads <- list(                  #Include these populations as reads
    "SHG"=c("SwedenSkoglund_MHG", "Motala_HG", "SwedenSkoglund_NHG"),
    "WHG"=c("SpanishMesolithic", "HungaryGamba_HG"), #SpanishMesolithic is the high coverage LaBrana1 I0585 and LaBrana2
    "EN"=c("LBK_EN", "HungaryGamba_EN", "Spain_EN", "Starcevo_EN", "LBKT_EN"),
    "MN"=c( "Spain_MN", "Baalberge_MN", "Iceman", "Esperstedt_MN" ),
    ## "Yamnaya"=c("Yamnaya"),
    ## "Poltavka"=c("Poltavka"),
    ## "Srubnaya"=c("Srubnaya"),
    "SteppeBA"=c("Yamnaya", "Poltavka", "Srubnaya"),
    "LN/BA"=c( "HungaryGamba_CA", "Alberstedt_LN", "Corded_Ware_LN", "Bell_Beaker_LN", "BenzigerodeHeimburg_LN", "Unetice_EBA", "HungaryGamba_BA", "Halberstadt_LBA")
    )
include.counts <- list(                 #Include these populations as hard calls. 
    "WHG"="Loschbour",
    "EN"="Stuttgart"
     )

mod.pops <- c("CEU", "GBR", "IBS", "TSI")
}
########################################################################
if(version=="v8"){
    map <- read.table(paste0("~/selection/code/files/v8/selection_label_", map ,"_map.txt"), as.is=TRUE, header=FALSE)
    int.include <- map[,2]
    names(int.include) <- map[,1]

    if("WHG" %in% int.include){
        int.names <- c("WHG", "SHG", "AEN", "CEM", "INC", "CLB", "STP")
        cols <- c("#4DAF4A","#4DAF4A", "#377EB8", "#377EB8", "#377EB8",  "#E41A1C",  "#E41A1C")
        include.counts <- list(                 #Include these populations as hard calls. 
                               "WHG"="Loschbour",
                               "CEM"="Stuttgart")
    }else{
        int.names <- c("HG", "AEN", "CEM", "INC", "CLB", "STP")
        cols <- c("#4DAF4A", "#377EB8", "#377EB8", "#377EB8",  "#E41A1C",  "#E41A1C")
        include.counts <- list(                 #Include these populations as hard calls. 
                               "HG"="Loschbour",
                               "CEM"="Stuttgart"
                               )
    }
    names(cols) <- int.names
    include.reads <- list()
    for(nm in int.names){include.reads[[nm]] <- names(int.include[int.include==nm])}

    
    mod.pops <- c("CEU", "GBR", "IBS", "TSI")
}

########################################################################

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- counts[,6:NCOL(counts)]
totals <- totals[,6:NCOL(totals)]

## Readme data and filter count data
readme <- read.table(readmefile, as.is=TRUE)
rownames(readme) <- readme[,2]
rownames(data) <- data[,1]
rownames(counts) <- data[,1]
rownames(totals) <- data[,1]
data <- data[rownames(readme),]
counts <- counts[rownames(readme),]
totals <- totals[rownames(readme),]

empty.data <- make.empty.data(int.names)
include.read.samples <- read.samples(indfile, include.reads)

########################################################################

n.plot.row <- floor((NROW(data)+1)/2)
n.plot.col <- ifelse(NROW(data)>1, 2,1)
pdf(paste0(out, outname), width=4*n.plot.col, height=3*n.plot.row+0.5*(n.plot.row==1))
par(mfrow=c(n.plot.row, n.plot.col))

mod.cols <- brewer.pal(4, "Set2")

for(i in 1:NROW(data)){
  plot(0,0, col="white", xlim=c(0.5, length(int.names)+0.5 ), ylim=ylim, bty="n", xlab="", ylab="Derived allele frequency", main=paste0(readme[i,1], " (", data$ID[i], ")"),xaxt="n", cex.axis=1.2, cex.main=1.2)
  cat(paste0("\r", readme[i,2]))
  chr <- data[i,"CHR"]
  read.data <- read.table(paste0(read.root, "jj2.chr", chr, ".readcounts.gz"), as.is=TRUE)
  read.data <- read.data[read.data[,1]==data[i,"ID"],]
    
  freq.data <- make.freq.data(int.names, include.reads, include.read.samples, include.counts,
                                read.data, counts[i,], totals[i,], empty.data)
  f <- fit.unconstrained.model.reads(freq.data, error.prob=error.prob)$p
  names(f) <- int.names
  eff.totals <- round(effective.data.size(freq.data),1)

  ## mod.f=range(counts[i,mod.pops]/totals[i,mod.pops])
  mod.f=counts[i,mod.pops]/totals[i,mod.pops]

  ci <- ci.unconstrained.model.reads(freq.data, error.prob, alpha=0.95)
  uci=ci$uci
  lci=ci$lci

  ## Flip to derived allele. 
  if(readme[data[i,1],3]==data[i,4]){
    f=1-f
    mod.f=1-mod.f
    uci=1-uci
    lci=1-lci
  }

  ## Plot
  segments(1:length(lci), lci, 1:length(uci), uci, lwd=4, col=cols)
  points(f, pch=16, cex=2, col=cols)
  text(f, labels=eff.totals, pos=2, cex=1.4)
  for(k in 1:length(mod.f)){
    abline(h=mod.f[k], lty=2, lwd=2, col=mod.cols[k])
    text( length(int.names)+0.5, mod.f[k], mod.pops[k], cex=1.4)
  }

  mtext(int.names, side=1, cex=1, line=1, at=1:length(int.names), las=2)
}
  
dev.off()

