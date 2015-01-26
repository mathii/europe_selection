#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
source("~/selection/code/lib/3pop_lib.R")
dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

#Modern GBR, CEU, FIN, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
root <- "~/selection/series/counts/all_3pop"
snproot <- "~/selection/series/snps/all_3pop"
## root <- "~/selection/series/counts/diet"
## snproot <- "~/selection/series/snps/diet"
out <- "~/selection/series/series/"
popfile="~/data/v6/hov6_1kg_motalav3_3pop.ind"

lambda=NA

########################################################################

gss <- c(100, 200)                               #Generations of selection
ss <- 10^(seq(log10(0.002), log10(0.1), length.out=10)) #Selection coefficient
Ne <- 6000                                          #2 Population size
N <- 1000                                           #Number of replicates
sig <- 10^-6.5                                      #genome-wide significance level

########################################################################


#Compute the likelihood 

## pops <- c("WHG", "Eneo", "YamnayaEBA", "CEU", "GBR", "IBS", "TSI")
## A <- matrix(c(0.187, 0.312, 0.501, 0.160, 0.413, 0.427, 0, 0.764, 0.236, 0, 0.714, 0.286),3, 4)

## Excluding IBS which seem to be messed up. 
pops <- c("WHG", "Eneo", "YamnayaEBA", "CEU", "GBR", "TSI")
A <- matrix(c(0.187, 0.312, 0.501, 0.160, 0.413, 0.427, 0, 0.714, 0.286),3, 3)

#Adding other ancient pops. 
## pops <- c("WHG", "Eneo", "YamnayaEBA", "CEU", "GBR", "TSI", "Spain_MN", "UneticeEBA", "GermanBellBeakerNeolithic", "CordedWareNeolithic")
## A <- matrix(c(0.187, 0.312, 0.501, 0.160, 0.413, 0.427, 0, 0.714, 0.286,
##               0.175, 0.825, 0.000, 0.328, 0.208, 0.464, 0.184, 0.365, 0.451,
##               0.243, 0.078, 0.679),nrow=3)


counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- data.matrix(counts[,6:NCOL(counts)])
totals <- data.matrix(totals[,6:NCOL(totals)])
counts <- counts[,pops]
totals <- totals[,pops]

#Select sites with f0 < 0.1
tc <- matrix(0, nrow=N, ncol=NCOL(counts))
tt <- matrix(0, nrow=N, ncol=NCOL(totals))
i=1
while(i <= N){
    try <- sample(NROW(counts),1)
    fr <- counts[try,4:NCOL(counts)]/totals[try,4:NCOL(totals)]
    fr <- fr[!is.infinite(fr)]
    if(mean(fr)>0.1){next}

    tc[i,] <- counts[try,]
    tt[i,] <- totals[try,]
    i <- i+1
}

results.all <- matrix(0,nrow=length(ss), ncol=length(gss))
results.one <- matrix(0,nrow=length(ss), ncol=length(gss))
#Add selection (to all pops)
for( gsi in 1:length(gss)){
    gs <- gss[gsi]
    for(si in 1:length(ss)){
        cat(paste0(gss[gsi]," ",  ss[si], "\n"))

        #All pops
        this.tc <- tc
        this.tt <- tt
        for(i in 1:N){
            for(j in 4:NCOL(counts)){
                this.fr <- pmax(0.01, this.tc[i,j]/this.tt[i,j])
                traj <- simulate.wright.fisher(Ne, gss[gsi], this.fr, ss[si])
                new.fr <- rev(traj)[1]
                this.tc[i,j] <- rbinom(1, this.tt[i,j], new.fr)
            }
            test <- test.3pop(this.tt[i,], this.tc[i,], A)
            if(test[2]<sig){results.all[si,gsi] <- results.all[si,gsi]+1} 
        }

        #One pop
        this.tc <- tc
        this.tt <- tt
        for(i in 1:N){
            j <- sample(4:NCOL(counts), 1)
            this.fr <- pmax(0.01, this.tc[i,j]/this.tt[i,j])
            traj <- simulate.wright.fisher(Ne, gs, this.fr, ss[si])
            new.fr <- rev(traj)[1]
            this.tc[i,j] <- rbinom(1, this.tt[i,j], new.fr)
            test <- test.3pop(this.tt[i,], this.tc[i,], A)
            if(test[2]<sig){results.one[si,gsi] <- results.one[si,gsi]+1} 
        }

    }
}
results.all <- results.all/N
results.one <- results.one/N


plot(ss, results.all[,1], col="#377EBA", type="b", pch=16, bty="n", lty=2, ylim=c(0,1), log="x", xlab="Selection coefficient", ylab="power")
lines(ss, results.all[,2], col="#E41A1C", type="b", pch=16, lty=2)
lines(ss, results.one[,1], col="#377EBA", type="b", pch=1, lty=3)
lines(ss, results.one[,2], col="#E41A1C", type="b", pch=1, lty=3)
legend("topleft", c("Selected in all populations", "Selected in one population", "100 generations of selection", "200 generations of selection"), col=c("black", "black", "#377EBA", "#E41A1C"), pch=16, lty=c(2,3,1,1), bty="n")
