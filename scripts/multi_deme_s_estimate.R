## Fit a multi deme model to a selection trajectory and
## estimate selection coefifficents
source("~/selection/code/lib/readlib.R")
dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

######################################################################

root <- "~/selection/counts/all"
## demes <- c("Spain", "Germany", "Hungary")
demes <- c("Spain", "Germany")
polynames.file <- "~/data/v6/use/polymap.txt"
modern <- c("IBS"="Spain", "CEU"="Germany")
datefile <- "~/data/v6/use/population_dates.txt"
gen.time <- 29
Ne <- 3000                             #TODO!
m <- 0.005
fix.m <- TRUE

snp <- "rs12913832"
flip <- TRUE

######################################################################

polynames <- read.table(polynames.file, as.is=TRUE)
population.demes <- modern
for(i in 1:NROW(polynames)){
    bits <- strsplit(polynames[i,2], "_")
    cnt <- bits[[1]][1]
    if(cnt %in% demes){
        population.demes <- c(population.demes, cnt)
        names(population.demes)[length(population.demes)] <-  polynames[i,1]
    }
}

## Majority call info
rd <- read.counts.and.data(root)
counts <- rd$counts
totals <- rd$totals
data <- rd$data

## Dates
dates <- read.table(datefile, as.is=TRUE, header=FALSE)
dd <- dates[,2]*1000
names(dd) <- dates[,1]
dd <- dd[names(dd) %in% names(population.demes)]
max.gen.ago <- ceiling(max(dd)/gen.time)


## This SNP
include <- data$ID==snp
snp.counts <- counts[include]
snp.totals <- totals[include]
names(snp.counts) <- names(snp.totals) <- colnames(counts)

obs <- list(N=array(0,c(1,length(demes),max.gen.ago)), N.A=array(0,c(1,length(demes),max.gen.ago)))
for(i in 1:length(population.demes)){
    what <- names(population.demes)[i]
    gen <- max.gen.ago-round(dd[what]/gen.time)
    where <- which(demes==population.demes[i])
    obs$N[1,where,gen] <- snp.totals[what]
    if(flip){
        obs$N.A[1,where,gen] <- snp.totals[what]-snp.counts[what]    
    }else{
        obs$N.A[1,where,gen] <- snp.counts[what]    
    }
}


if(fix.m){
    est<-estimate.s.m(obs, Ne, M=m, update="Soft EM", max.iters=10, verbose=TRUE)
    est.fix<-estimate.s.m(obs, Ne, M=m, update="Soft EM", max.iters=10, verbose=TRUE, s.free=FALSE)
}else{
    est<-estimate.s.m(obs, Ne, M=NULL, update="Soft EM", max.iters=10, verbose=TRUE, initial.M=m)
    est.fix<-estimate.s.m(obs, Ne, M=NULL, update="Soft EM", max.iters=10, verbose=TRUE, initial.M=m, s.free=FALSE)
}
par(mfrow=c(1,2))
plot.wright.fisher.lattice.observations(obs, est$f, est$f, est.s=est$s, error.bars=TRUE, main="s estimate free")
plot.wright.fisher.lattice.observations(obs, est.fix$f, est.fix$f, est.s=est.fix$s, error.bars=TRUE, main="s estimate fixed")

X2 <- 2*(est$log.likelihood-est.fix$log.likelihood)
p <- pchisq(X2, df=length(demes)-1, lower.tail=FALSE)

print(est$s)
print(est$M)
print(est$log.likelihood)

print(est.fix$s)
print(est.fix$M)
print(est.fix$log.likelihood)

## free <- rep(0,10)
## fixed <- rep(0,10)
## ms <- seq(0.001,0.01,0.001)
## free.store <- list()
## fixed.store <- list()
## for(i in 1:10){
##     est<-estimate.s.m(obs, Ne, M=ms[i], update="Soft EM", max.iters=10, verbose=TRUE)
##     est.fix<-estimate.s.m(obs, Ne, M=ms[i], update="Soft EM", max.iters=10, verbose=TRUE, s.free=FALSE)
##     free[i] <- est$log.likelihood
##     fixed[i] <- est.fix$log.likelihood
##     free.store[[i]] <- est
##     fixed.store[[i]] <- est.fix
## }
