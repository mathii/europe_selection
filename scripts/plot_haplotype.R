#Make a haplotype plot around a SNP
source("~/selection/code/lib/readlib.R")

############################################################
snp <- "rs3827760"
flank <- 100000
root <- "~/selection/counts/all"
pops <- c("CHB", "Motala_HG", "CEU")
out <- "~/selection/analysis/series/EDAR_haplotype.pdf"
read.root <- "~/data/v6/reads/jj2"

pops.with.reads <- "Motala_HG"
cols <- c("blue", "purple", "red")

############################################################

## Load all snp info and select SNPs to be used 
rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

this.snpinfo <- data[data[,1]==snp,]
include <- (data[,2]==this.snpinfo[,2])&(abs(data[,3]-this.snpinfo[,3])<flank)
data <- data[include,]
counts <- counts[include,]
totals <- totals[include,]


