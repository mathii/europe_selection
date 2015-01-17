#Estimate selection coefficients

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
include.ci <- TRUE

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

results <- matrix(NA, ncol=2, nrow=NROW(data))
if(include.ci){results <- matrix(NA, ncol=4, nrow=NROW(data))}
rownames(results) <- data$ID

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
            this.est <- estimate.s(this.traj, Ne, method="Soft EM", verbose=FALSE)
            this.p <- p.value(this.traj, Ne, this.est$s)
            if(include.ci){
                ci <- find.confidence.interval(this.traj, Ne, this.est$s)
                results[i,] <- c(this.est$s, this.p, ci)
            } else{
                results[i,] <- c(this.est$s, this.p)        
            }
        }
        , error=function(e){e}
        )
}

write.table(results, paste0("~/selection/analysis/s_estimates/s_estimates_", tag, ".txt"), row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
