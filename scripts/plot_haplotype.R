#Make a haplotype plot around a SNP
source("~/selection/code/lib/readlib.R")

############################################################
snp <- "rs3827760"
flank <- 100000
root <- "~/data/v6/use/v61kg_europe2names"
pops <- c("CHB", "Motala_HG", "CEU")
out <- "~/selection/analysis/EDAR/EDAR_haplotype"
read.root <- "~/data/v6/reads/jj2"

pops.with.reads <- "Motala_HG"
cols <- c("blue", "purple", "red")

############################################################

## Load all snp info and select SNPs to be used 
data <- read.table(paste0(root, ".snp"))
this.snpinfo <- data[data[,1]==snp,]
include <- (data[,2]==this.snpinfo[,2])&(abs(data[,4]-this.snpinfo[,4])<flank)
data <- data[include,]

write.table(data, paste0(out, ".snp"), col.names=FALSE, row.names=FALSE, quote=FALSE)

## At this point, stop and run the following command to pull out the CHB haploypes:
## python ~/spindrift/Freq.py -d ~/data/v6/use/v61kg_europe2names -p CHB,CEU -o CHB_CEU_EDAR_haps -s EDAR_haplotype -g

ind <- read.table(paste0(root, ".ind"), as.is=TRUE)
gt <- read.table("~/selection/analysis/EDAR/CHB_CEU_EDAR_haps.gt", as.is=TRUE, header=TRUE)
gt.data <- gt[,1:5]
gt <- gt[,6:NCOL(gt)]

## Load reads
reads <- read.table(paste0(read.root, ".chr", this.snpinfo[,2], ".readcounts"), as.is=TRUE, header=FALSE)
reads <- reads[reads[,1] %in% gt.data[,1],]

ntotal <- 0
for(p in pops){
    ntotal <- ntotal+sum(ind[,3]==p)
}

ht <- matrix(0, nrow=ntotal, ncol=NROW(gt))
pr <- matrix(0, nrow=ntotal, ncol=NROW(gt))

i=1
for(pop in pops){
    if(pop %in% pops.with.reads){
        for(indiv in ind[ind[,3]==pop,1]){
            j=1
            for(snpid in gt.data[,1]){
                select <- reads[,1]==snpid & reads[,2]==indiv
                ref=reads[select,3]
                alt=reads[select,4]
                if(alt==0 & ref==0){    #No data
                    ht[i,j] <- 1
                    pr[i,j] <- 0
                }else if(alt==0){
                    ht[i,j] <- 2
                    pr[i,j] <- 1-0.5^(ref-1)
                }else if(ref==0){
                    ht[i,j] <- 0
                    pr[i,j] <- 1-0.5^(alt-1)
                }else{                  #Het
                    ht[i,j] <- 1
                    pr[i,j] <- 1
                }
                j=j+1
            }
            i=i+1
        }
    }else{
        for(indiv in ind[ind[,3]==pop,1]){
            ht[i,] <- gt[,indiv]
            pr[i,] <- 1
            i=i+1
        }
    }
}
