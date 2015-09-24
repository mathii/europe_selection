#Make a haplotype plot around a SNP
source("~/selection/code/lib/readlib.R")
library(RColorBrewer)

############################################################
## snp <- "rs3827760"
## flank <- 150000
## range <- NA
## what<-"EDAR"
## pops <- c("CHB", "Motala_HG", "CEU")
## pops.with.reads <- "Motala_HG"
## mod.pops <- "CHB<CEU"

############################################################
## snp <- "rs1426654"
## flank <- 150000
## range <- NA
## what<-"SLC24A5"
## pops <- c("CHB", "Motala_HG", "CEU")
## pops.with.reads <- "Motala_HG"
## version <- "v8"
## mod.pops <- c("CHB","CEU")

############################################################
## snp <- "rs4988235"
## flank <- 150000
## range <- NA
## what <- "LCT"
## pops.with.reads <- c("Motala_HG","Yamnaya_Samara", "Sweden_PWC", "Srubnaya", "Poltavka", "EHG", "Srubnaya_1d_rel_I0430", "Central_LNBA", "Bell_Beaker_LN", "Yamnaya_Kalmykia","Northern_LNBA", "Hungary_BA")
## pops <- c("CEU", pops.with.reads)
## mod.pops <- "CEU"

############################################################

snp <- NA
flank <- 150000
range <- c(2,46524541,46613842)
what<-"EPAS1"
pops.with.reads <- c("Botigiriayocc_LIP","Huaca_Prieta_IP","Inca_LH","Lauricocha_EMA","Lauricocha_IP","Lauricocha_LA","Lima_EIP","Nasca_EIP","Pacapaccari_LIP","Tiwanaku_MH","Wari_Coast_MH","Wari_Highlands_MH","Ychsma_LIP")
version <- "peru"
mod.pops <- c("CHB","PEL")
pops <- c(mod.pops[1], pops.with.reads, mod.pops[2])

############################################################

root <- paste0("~/data/",version,"/use/",version,"1kg_europe2names")
read.root <- paste0("~/data/",version,"/reads/jj2")
subsample <- TRUE
cols <- brewer.pal(3, "Set1")

############################################################

out <- paste0("~/selection/analysis/",version,"/", what,"/",what ,"_haplotype")

## Load all snp info and select SNPs to be used 
data <- read.table(paste0(root, ".snp"))
if(all(is.na(range))){
    this.snpinfo <- data[data[,1]==snp,]
    include <- (data[,2]==this.snpinfo[,2])&(abs(data[,4]-this.snpinfo[,4])<flank)
    this.chr <- this.snpinfo[,2]
}else{
    include <- data[,2]==range[1] & data[,4]>range[2] & data[,4]<range[3]
    this.chr <- range[1]
}
data <- data[include,]

write.table(data, paste0(out, ".snp"), col.names=FALSE, row.names=FALSE, quote=FALSE)

## At this point, stop and run the following command to pull out the CHB haploypes:
## python ~/spindrift/Freq.py -d ~/data/${V}/use/${V}1kg_europe2names -p ${mod.pops} -o test_haps -s ${what}_haplotype -g
if(!file.exists(paste0("~/selection/analysis/",version,"/",what,"/test_haps.gt"))){
    stop("At this point...")
}

ind <- read.table(paste0(root, ".ind"), as.is=TRUE)
gt <- read.table(paste0("~/selection/analysis/",version,"/",what,"/test_haps.gt"), as.is=TRUE, header=TRUE)

gt.data <- gt[,1:5]
gt <- gt[,6:NCOL(gt)]

## Load reads
reads <- read.table(paste0(read.root, ".chr", this.chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
reads <- reads[reads[,1] %in% gt.data[,1],]

ntotal <- 0
for(p in pops){
    ntotal <- ntotal+sum(ind[,3]==p)
}

ht <- matrix(0, nrow=ntotal, ncol=NROW(gt))
pr <- matrix(0, nrow=ntotal, ncol=NROW(gt))

snpi=1
poplist=c()
indivlist=c()
for(snpid in gt.data[,1]){
    cat(paste0("\r", snpi, "/", NROW(gt.data)))
    this.reads <- reads[reads[,1]==snpid,]
    indivi=1
    for(pop in pops){
        if(pop %in% pops.with.reads){
            for(indiv in ind[ind[,3]==pop,1]){
                if(snpi==1){
                    poplist <- c(poplist, pop)
                    indivlist <- c(indivlist, indiv)
                }
                select <-  this.reads[,2]==indiv
                ref=this.reads[select,3]
                alt=this.reads[select,4]
                if(alt==0 & ref==0){    #No data
                    ht[indivi,snpi] <- NA        #Using NA to mean missing. 
                    pr[indivi,snpi] <- 0
                }else if(alt==0){
                    ht[indivi,snpi] <- 2
                    pr[indivi,snpi] <- 1-0.5^(ref-1)
                }else if(ref==0){
                    ht[indivi,snpi] <- 0
                    pr[indivi,snpi] <- 1-0.5^(alt-1)
                }else{                  #Het
                    ht[indivi,snpi] <- 1
                    pr[indivi,snpi] <- 1
                }
                indivi=indivi+1
            }
        } else{
            for(indiv in ind[ind[,3]==pop,1]){
                if(snpi==1){
                    poplist <- c(poplist, pop)
                    indivlist <- c(indivlist, indiv)
                }
                ht[indivi,snpi] <- gt[snpi,indiv]
                pr[indivi,snpi] <- 1
                indivi=indivi+1
            }
        }
    }
    snpi=snpi+1
}

## Flip to European alleles
if(what=="EDAR"){
    for(i in 1:NCOL(ht)){
        if(mean(ht[poplist=="CEU",i], na.rm=TRUE)>1){ht[,i] <- 2-ht[,i]}
    }
}

if(what=="EDAR"){
    plot.order <- c("I0011", "I0012", "I0017", "I0013", "I0014", "I0015")
} else if(what=="LCT"){
    plot.order <- rev(c("I0232","I0011","I0357","Ajvide52","I0424","I0444","I0126","I0211","I0235","I0361","I0421","I0431","I0441","I0804","I1544","I1546","RISE240","RISE431","RISE435","RISE546","RISE550","RISE98","I1504","I0164","I0112","I0430","I0423"))
} else{
    plot.order <- indivlist[poplist %in% pops.with.reads]
}
    ## No subsampling
## subsample 20 CEU and 20 CHB
if(subsample){
    if(length(mod.pops)==1){
        s.size <- 40
        sub <- c(which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], sample(which(poplist==mod.pops[1]), s.size))
    }else if(length(mod.pops)==2){
        s.size <- 20
        sub <- c(sample(which(poplist==mod.pops[1]), s.size), which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], sample(which(poplist==mod.pops[2]), s.size))
    }else{
        stop("px1")
    }
}else{
    if(length(mod.pops)==1){
        s.size <- sum(poplist==mod.pops[1])
        sub <- c(which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], sample(which(poplist==mod.pops[1]), s.size))
    }else if(lengt(mod.pops)==2){
        s.size.1 <- sum(poplist==mod.pops[1])
        s.size.2 <- sum(poplist==mod.pops[2])
        sub <- c(sample(which(poplist==mod.pops[1]), s.size.1), which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], sample(which(poplist==mod.pops[2]), s.size.2))
    } else{
        stop("px2")
    }
}
sub.ht <- ht[sub,]
sub.pr <- pr[sub,]
subtotal <- length(sub)

remove <- rep(FALSE, NCOL(sub.ht))
## for(i in 1:NCOL(sub.ht)){
##     remove[i] <- all(sub.ht[sub%in%which(poplist=="CEU"),i]==0)|all(sub.ht[sub%in%which(poplist=="CEU"),i]==2)
## }
## sub.ht <- 2-sub.ht

sub.ht <- sub.ht[,!remove]
sub.pr <- sub.pr[,!remove]
if(!(is.na(snp))) {
    snpi <- which(gt.data[!remove,1]==snp)
} else{
    snpi <- NA
}

## cols <- c("grey", "pink", "darkred", "white")
## s.cols <- c("grey", "lightblue", "darkblue", "white")

cols <- c("grey", "pink", "darkred", "white")
s.cols <- c("grey", "lightblue", "darkblue", "white")

nsnp <- NROW(gt)

if(subsample){
    pdf(paste0(out, "_subsampled.pdf"), width=round(nsnp/10), height=round(subtotal/10))
}else{
    pdf(paste0(out, ".pdf"), width=12, height=24)
}
plot(0,0, col="white", bty="n", xaxt="n", yaxt="n", xlim=c(0,nsnp), ylim=c(0,subtotal), xlab="", ylab="")

   for(i in 1:subtotal){
    cc <- cols[sub.ht[i,]+1]
    bc=col2rgb(cols[sub.ht[i,]+1], alpha=TRUE)
    bc[4,] <- round(255*sub.pr[i,]^2)
    bc <- apply(bc, 2, function(x)do.call(rgb, as.list((x/255))))
    points(1:nsnp, rep(i,nsnp), pch=21, cex=1, col=cc, bg=bc)
}

if(!is.na(snpi)){
cc=s.cols[sub.ht[,snpi]+1]
bc=col2rgb(s.cols[sub.ht[,snpi]+1], alpha=TRUE)
bc[4,] <- round(255*sub.pr[,snpi]^2)
bc <- apply(bc, 2, function(x)do.call(rgb, as.list((x/255))))
points(rep(snpi,subtotal), 1:subtotal, pch=21, cex=1, col=cc, bg=bc)
}

if(what=="LCT"){
    abline(h=subtotal-s.size+0.5, col="black")
    abline(h=6.5, col="black")
    abline(h=1.5, col="black")
    abline(h=0.5, col="black")
    abline(h=21.5, col="black")
}
if(what!="LCT"){
    abline(h=s.size+sum(!(poplist %in% mod.pops))+0.5, col="black")
    abline(h=s.size+0.5, col="black")
    mtext(mod.pops, side=2, at=c(s.size/2, 1.5*s.size+sum(!(poplist %in% mod.pops))), line=-1, las=2)
    mtext("Ancient", side=2, at=c( s.size+0.5*sum(!(poplist %in% mod.pops))), line=-1, las=2)

}
dev.off()

