#Test whether the test is robust to changing the admixture proportions.
source("~/selection/code/lib/3pop_lib.R")
source("~/Packages/s_lattice/simulation.R")

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

########################################################################

rnds <- seq(0,1,length.out=11)                               #proportions of randomness
Nsims <- NA                                           #Number of replicates per test
Npowersims <- 10000                                           #Number of replicates per test
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


degf <- dim(A)[2]

counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- data.matrix(counts[,6:NCOL(counts)])
totals <- data.matrix(totals[,6:NCOL(totals)])
counts <- counts[,pops]
totals <- totals[,pops]


#Select sites with f0 < 0.1
tc <- counts
tt <- totals
if(!is.na(Nsims)){
    s.witch <- sample(NROW(counts), Nsims)
    counts <- counts[s.witch,]
    totals <- totals[s.witch,]
}

lambda.all <- rep(0, length(rnds))
#Add selection (to all pops)
for( rndi in 1:length(rnds)){
    rnd <- rnds[rndi]
    cat(paste0(rnd, "\n"))
    this.res <- matrix(0, nrow=NROW(totals), ncol=2)
    
    for(i in 1:NROW(totals)){
        N <- totals[i,]
        N.A <- counts[i,]

        if(all(N==N.A)){
            this.res[i,] <- NA
            next
        }

        Atmp <- A
        for(j in 1:NCOL(A)){
            kg <- rexp(3)
            rndp <- kg/sum(kg)
            Atmp[,j] <- rnd*rndp+(1-rnd)*A[,j]
        }
        
        this.res[i,] <- test.3pop(N, N.A, Atmp)
}
    lambda.all[rndi] <- median(this.res[!is.na(this.res[,2]),1])/qchisq(0.5, df=degf)
}

## pdf("~/selection/paper/FigureS4a.pdf")
plot(rnds, lambda.all, col="#377EBA", type="b", pch=16, bty="n", lwd=1, xlab="Random proportion", ylab="Genomic inflation factor")
## dev.off()

#Now test power.
gens <- 200
s <- 0.01
Ne <- 6000

tc <- matrix(0, nrow=Npowersims, ncol=NCOL(counts))
tt <- matrix(0, nrow=Npowersims, ncol=NCOL(totals))
i=1
while(i <= Npowersims){
    try <- sample(NROW(counts),1)
    fr <- counts[try,4:NCOL(counts)]/totals[try,4:NCOL(totals)]
    fr <- fr[!is.infinite(fr)]
    if(mean(fr)>0.1){next}

    tc[i,] <- counts[try,]
    tt[i,] <- totals[try,]
    i <- i+1
}

one.power <- c(0, length(rnds))
for( rndi in 1:length(rnds)){
    rnd <- rnds[rndi]
    cat(paste0(rnd, "\n"))

    results.one <- 0
    
    this.tc <- tc
    this.tt <- tt
    for(i in 1:Npowersims){
        ## Selection in one population
        j <- sample(4:NCOL(counts), 1)
        this.fr <- pmax(0.01, this.tc[i,j]/this.tt[i,j])
        traj <- simulate.wright.fisher(Ne, gens, this.fr, s)
        new.fr <- rev(traj)[1]
        this.tc[i,j] <- rbinom(1, this.tt[i,j], new.fr)

        ## Randomise A
        Atmp <- A
        for(j in 1:NCOL(A)){
            kg <- rexp(3)
            rndp <- kg/sum(kg)
            Atmp[,j] <- rnd*rndp+(1-rnd)*A[,j]
        }
        
        test <- test.3pop(this.tt[i,], this.tc[i,], Atmp)
        corrected.p <- pchisq(test[1]/lambda.all[rndi], df=degf, lower.tail=F)
        if(corrected.p<sig){results.one <- results.one+1} 
    }
    one.power[rndi] <- results.one/Npowersims


    

}
dev.new()
## pdf("~/selection/paper/FigureS4b.pdf")
plot(rnds, one.power, col="#E41A1C", type="b", pch=16, bty="n", lty=3, xlab="Random proportion", ylab="Power")
## dev.off()
