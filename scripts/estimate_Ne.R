###############################################################

root <- "~/selection/series/counts/random_sample_1k_spanish_german_noWHG_v6"
datefile <- "~/selection/data/simple_dates.txt"
max.gen.ago <- 250
gen.time <- 29

###############################################################

dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

###############################################################

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- counts[,6:NCOL(counts)]
totals <- totals[,6:NCOL(totals)]

dates <- read.table(datefile, as.is=TRUE, header=TRUE)
dates$mid <- (dates[,5]+dates[,6])/2
dates <- dates[,c("Pop", "mid")]
dates<-aggregate(dates[,2], by=list(dates[,1]), FUN=mean)
dd <- dates[,2]
names(dd) <- dates[,1]
df <- rep(0,5)
names(df) <- c("CEU", "GBR", "FIN", "IBS", "TSI")
dates <- c(dd, df)

log10Negrid <- seq(3,4,length.out=20)

results <- matrix(0, ncol=length(log10Negrid), nrow=NROW(data))
rownames(results) <- data$ID

for(i in 1:NROW(data)){
    cat(paste0("\r", i))
    this.traj <- data.frame(N=rep(0,max.gen.ago), N.A=rep(0,max.gen.ago))
    for(j in 1:NCOL(counts)){
        what <- names(counts)[j]
        gen <- max.gen.ago-round(dates[what]/29)
        if(is.na(gen)){stop(paste0("No date for ", what))}
        this.traj$N[gen] <- this.traj$N[gen]+totals[i,j]
        this.traj$N.A[gen] <- this.traj$N.A[gen]+counts[i,j]
    }

    for(j in 1:length(log10Negrid)){
        Ne <- 10^log10Negrid[j]
        results[i,j] <- wfhmm.call(this.traj, 0, Ne, viterbi=FALSE, paths=0, likelihood="Model")$log.likelihood
    }
}

pdf("~/selection/paper/Figure_SA.pdf")
profile <- colSums(results-results[,1])
fitspln <- smooth.spline(10^log10Negrid, profile)
xpts <- seq(4000,10000, length.out=100)
pts <- predict(fitspln, xpts)
plot(pts, log="x", type="l", col="blue", ylab="log-likelihood (+constant)", xlab=expression('2N'[e]))

drv.fun <- approxfun(predict(fitspln, xpts, deriv=1))
ll.max <- uniroot(drv.fun, interval=c(4000, 10000))$root
print(ll.max)
abline(v=ll.max, col="red") 

act.fun <- approxfun(predict(fitspln, xpts))
act.fun <- approxfun(xpts, predict(fitspln, xpts)$y-act.fun(ll.max)+qchisq(0.95, df=1)/2)
lo<-uniroot(act.fun, c(4000, ll.max))
hi<-uniroot(act.fun, c(ll.max, 10000))
print(c(lo, hi))
abline(v=lo, col="red", lty=2) 
abline(v=hi, col="red", lty=2) 
dev.off()
