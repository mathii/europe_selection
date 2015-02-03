#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
source("~/selection/code/lib/3pop_lib.R")
library(reshape2)
dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

#Modern GBR, CEU, FIN, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## Details
root <- "~/selection/counts/all"
snproot <- "~/data/v6/use/v61kg_europe2names"
out <- "~/selection/analysis/power/"
indfile <- "~/data/v6/use/v61kg_europe2names.ind"

lambda=1.20

########################################################################

gss <- c(50, 100, 200)                               #Generations of selection
ss <- 10^(seq(log10(0.002), log10(0.1), length.out=10)) #Selection coefficient
Ne <- 6000                                          #2 Population size
N <- 1000                                           #Number of replicates
sig <- 10^-6.79                                      #genome-wide significance level

pops <- c("WHG", "EN", "Yamnaya", "CEU", "GBR", "IBS", "TSI")
#Check if the SNP is monomorphic in these populations. 
A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.227, 0, 0.712, 0.288),3, 4) 
degf <- dim(A)[2]

########################################################################


counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)

counts$WHG <- counts$Loschbour+counts$LaBrana1+counts$HungaryGamba_HG
totals$WHG <- totals$Loschbour+totals$LaBrana1+totals$HungaryGamba_HG
counts$EN <- counts$LBK_EN +counts$Stuttgart+counts$HungaryGamba_EN+counts$Spain_EN+counts$Starcevo_EN+counts$LBKT_EN
totals$EN <- totals$LBK_EN +totals$Stuttgart+totals$HungaryGamba_EN+totals$Spain_EN+totals$Starcevo_EN+totals$LBKT_EN


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
            test[2] <- pchisq(test[1]/lambda, df=degf, lower.tail=FALSE)
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
            test[2] <- pchisq(test[1]/lambda, df=degf, lower.tail=FALSE)
            if(test[2]<sig){results.one[si,gsi] <- results.one[si,gsi]+1} 
        }

    }
}
results.all <- results.all/N
results.one <- results.one/N


pdf("~/selection/analysis/power/majority_power.pdf")
plot(ss, results.all[,1], col="#377EBA", type="b", pch=16, bty="n", lty=2, ylim=c(0,1), log="x", xlab="Selection coefficient", ylab="power")
lines(ss, results.all[,2], col="#E41A1C", type="b", pch=16, lty=2)
lines(ss, results.all[,3], col="#4DAF4A", type="b", pch=16, lty=2)
lines(ss, results.one[,1], col="#377EBA", type="b", pch=1, lty=3)
lines(ss, results.one[,2], col="#E41A1C", type="b", pch=1, lty=3)
lines(ss, results.one[,3], col="#4DAF4A", type="b", pch=16, lty=3)
legend("topleft", c("Selected in all populations", "Selected in one population", "50 generations of selection", "100 generations of selection", "200 generations of selection"), col=c("black", "black", "#377EBA", "#E41A1C", "#4DAF4A"), pch=16, lty=c(2,3,1,1,1), bty="n")
dev.off()

colnames(results.all) <- colnames(results.one) <- gss
m1 <- melt(results.all)
m1[,1] <- "All"
m2 <- melt(results.one)
m2[,1] <- "One"
mres <- rbind(m1, m2)
mres <- cbind(mres, "Majority")
colnames(mres) <- c("Populations", "Generations", "Power", "Analysis")
write.table(mres, paste0(out, "power.majority.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
