#Estimate log-likelihoods for the selection coefficient. 

args <- commandArgs(TRUE)
print(args)
from <- as.numeric(args[1])
to <- as.numeric(args[2])
tag <- args[3]

###############################################################

root <- "~/selection/counts/all"
datefile <- "~/data/v6/use/population_dates.txt"
include <- "~/data/v6/use/used_series.txt"
gen.time <- 29
Ne <- 6000                             #TODO!

###############################################################

s.points <- c(0, 10^seq(-4,-1, length.out=7))

###############################################################

dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

###############################################################

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)

counts <- counts[from:(min(NROW(counts),to)),]
totals <- totals[from:(min(NROW(totals),to)),]

data <- counts[,1:5]
counts <- counts[,6:NCOL(counts)]
totals <- totals[,6:NCOL(totals)]

include <- scan(include, what="")
counts <- counts[,include]
totals <- totals[,include]

dates <- read.table(datefile, as.is=TRUE, header=FALSE)
dd <- dates[,2]*1000
names(dd) <- dates[,1]
dd <- dd[include]

max.gen.ago <- ceiling(max(dd)/gen.time)

likelihoods <- matrix(NA, ncol=length(s.points), nrow=NROW(data))
rownames(likelihoods) <- data$ID

for(i in 1:NROW(data)){
    cat(paste0("\r", i))
    this.traj <- data.frame(N=rep(0,max.gen.ago), N.A=rep(0,max.gen.ago))
    for(j in 1:NCOL(counts)){
        what <- names(counts)[j]
        gen <- max.gen.ago-round(dd[what]/gen.time)
        if(is.na(gen)){stop(paste0("No date for ", what))}
        this.traj$N[gen] <- this.traj$N[gen]+totals[i,j]
        this.traj$N.A[gen] <- this.traj$N.A[gen]+counts[i,j]
    }
    #Ignore errors if something randomly goes wrong. 
    e <- tryCatch(
        if(sum(this.traj$N[1:(NROW(this.traj)-1)]>1)>2){
            these.likelihoods <- NA*s.points
            for(j in 1:length(s.points)){
                likelihoods[i,j] <- wfhmm.call(this.traj, s.points[j], Ne, h=0.5, viterbi=FALSE, paths=0, likelihood="Model")$log.likelihood
            }
        }
        , error=function(e){e}
        )
}

write.table(likelihoods, paste0("~/selection/analysis/s_estimates/s_likelihoods_", tag, ".txt"), row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
