## Test how affected the s_estimate is by admixture.

###############################################################

dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

###############################################################

g1 <- 100
g2 <- 200
ss <- 10^(seq(-3,-1,length.out=11))
as <- c(0,0.1,0.5)

N.Ne <- 100
N.s <- 100

Ne <- 10000
Ne.grid <- seq(1000,20000, length.out=10)
Ne.grid.fine <- seq(1000,20000, length.out=1000)

###############################################################
## Function: simulate a population that gets admixture from
## another pop. pop1 --a--> pop2

sim.freq.admix.pop <- function(Ne1, g1, a, Ne2=Ne1, g2=g1,s1=0,s2=s1,f01=runif(1),f02=f01){
    ## Up to the admixture
    traj1 <- simulate.wright.fisher(Ne1, g1, f01, s1)
    traj2 <- simulate.wright.fisher(Ne2, g1, f02, s2)
    ## admixture proportion
    mixf <- (1-a)*rev(traj1)[1]+a*rev(traj2)[1]
    ## rest of the trajectory
    traj3 <- simulate.wright.fisher(Ne2, g2-g1, mixf, s2)
    return(c(traj2, traj3))
}

###############################################################
## Function: simulate observations
## Two populations of size Ne, admixture a after g1 generations
## g2 in total. selection s in the receiving population.

sim.freq.admix.obs <- function(Ne, g1, g2, f0, a, s, sample.pattern){
    if(f0<=0|f0>=1){stop("f0 must be intermediate")}
    f <- 0
    while(rev(f)[1]<=0 | rev(f)[1]>=1){
        f <- sim.freq.admix.pop(Ne1=Ne, g1=g1, a=a, Ne2=Ne, g2=g2, s1=0, s2=s, f01=f0, f02=f0)
    }
    obs <- generate.observations.from.path(f, sample.pattern)
    return(obs)
}


###############################################################
## 1. Get MLE estimates of Ne, given alpha

Ne.MLE <- 0*as
Ne.ss <- list()
Ne.ll.fine <- list()
sample.pattern<-rep(0,200)
sample.pattern[c(20,40,60,80,100,120,140,160,180)]<-20
sample.pattern[200]<-200

for(i in 1:length(as)){
    a=as[i]
    ll <- matrix(0, nrow=N.Ne, ncol=length(Ne.grid))
    for(j in 1:N.Ne){
        cat(paste0("\r", a, ":",j, "/", N.Ne))
        obs <- sim.freq.admix.obs(Ne, g1, g2, runif(1), a, 0, sample.pattern)
        for(k in 1:length(Ne.grid)){
            ll[j,k] <- wfhmm.call(obs, 0, Ne.grid[k], viterbi=FALSE, paths=0, likelihood="Model")$log.likelihood
        }
    }
    ll <- colSums(ll)
    ss <- smooth.spline(Ne.grid, ll)
    Ne.ss[[i]] <- ss
    ll.fine <- predict(ss, Ne.grid.fine)
    Ne.ll.fine[[i]] <- ll.fine
    Ne.MLE[i] <- ll.fine$x[which.max(ll.fine$y)]
}

###############################################################
## 1. Using the MLE estimates of Ne, estimate selection coefficients. 
