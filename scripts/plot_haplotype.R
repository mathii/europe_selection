#Make a haplotype plot around a SNP
source("~/selection/code/lib/readlib.R")
library(RColorBrewer)

############################################################
## snp <- "rs3827760"
## flank <- 150000
## what<-"EDAR"
## pops <- c("CHB", "Motala_HG", "CEU")
## pops.with.reads <- "Motala_HG"

############################################################
snp <- "rs1426654"
flank <- 150000
what<-"SLC24A5"
pops <- c("CHB", "Motala_HG", "CEU")
pops.with.reads <- "Motala_HG"

############################################################
## snp <- "rs4988235"
## flank <- 150000
## what <- "LCT"
## pops.with.reads <- c("Motala_HG","Yamnaya_Samara", "Sweden_PWC", "Srubnaya", "Poltavka", "EHG", "Srubnaya_1d_rel_I0430", "Central_LNBA", "Bell_Beaker_LN", "Yamnaya_Kalmykia","Northern_LNBA", "Hungary_BA")
## pops <- c("CEU", pops.with.reads)
############################################################

root <- "~/data/v8/use/v81kg_europe2names"
read.root <- "~/data/v8/reads/jj2"
subsample <- TRUE
cols <- brewer.pal(3, "Set1")

############################################################

out <- paste0("~/selection/analysis/v8/", what,"/",what ,"_haplotype")

## Load all snp info and select SNPs to be used 
data <- read.table(paste0(root, ".snp"))
this.snpinfo <- data[data[,1]==snp,]
include <- (data[,2]==this.snpinfo[,2])&(abs(data[,4]-this.snpinfo[,4])<flank)
data <- data[include,]

write.table(data, paste0(out, ".snp"), col.names=FALSE, row.names=FALSE, quote=FALSE)

## At this point, stop and run the following command to pull out the CHB haploypes:
## python ~/spindrift/Freq.py -d ~/data/v6/use/v61kg_europe2names -p CHB,CEU -o test_haps -s ${what}_haplotype -g

ind <- read.table(paste0(root, ".ind"), as.is=TRUE)
gt <- read.table(paste0("~/selection/analysis/v8/",what,"/test_haps.gt"), as.is=TRUE, header=TRUE)

gt.data <- gt[,1:5]
gt <- gt[,6:NCOL(gt)]

## Load reads
reads <- read.table(paste0(read.root, ".chr", this.snpinfo[,2], ".readcounts.gz"), as.is=TRUE, header=FALSE)
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
    if(what!="LCT"){
        s.size <- 20
        sub <- c(sample(which(poplist=="CHB"), s.size), which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], sample(which(poplist=="CEU"), s.size))
    }else{
        s.size <- 40
        sub <- c(which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], sample(which(poplist=="CEU"), s.size))
    }
}else{
    s.size <- sum(poplist=="CHB")
    sub <- c(which(poplist=="CHB"), which(poplist %in% pops.with.reads)[match(plot.order, indivlist[poplist %in% pops.with.reads])], which(poplist=="CEU"))
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
snpi <- which(gt.data[!remove,1]==snp)


## cols <- c("grey", "pink", "darkred", "white")
## s.cols <- c("grey", "lightblue", "darkblue", "white")

cols <- c("grey", "pink", "darkred", "white")
s.cols <- c("grey", "lightblue", "darkblue", "white")

nsnp <- NROW(gt)

if(subsample){
    pdf(paste0(out, "_subsampled.pdf"), width=12, height=6)
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
cc=s.cols[sub.ht[,snpi]+1]
bc=col2rgb(s.cols[sub.ht[,snpi]+1], alpha=TRUE)
bc[4,] <- round(255*sub.pr[,snpi]^2)
bc <- apply(bc, 2, function(x)do.call(rgb, as.list((x/255))))
points(rep(snpi,subtotal), 1:subtotal, pch=21, cex=1, col=cc, bg=bc)

if(what=="LCT"){
    abline(h=subtotal-s.size+0.5, col="black")
    abline(h=6.5, col="black")
    abline(h=1.5, col="black")
    abline(h=0.5, col="black")
    abline(h=21.5, col="black")
}
if(what!="LCT"){
    abline(h=s.size+sum(poplist=="Motala_HG")+0.5, col="black")
    abline(h=s.size+0.5, col="black")
    mtext(c("CHB", "Motala_HG", "CEU"), side=2, at=c(s.size/2, s.size+sum(poplist=="Motala_HG")/2, 1.5*s.size+sum(poplist=="Motala_HG")), line=-2, las=2)
}
dev.off()

