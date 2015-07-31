#Test whether the test is robust to changing the admixture proportions.
source("~/selection/code/lib/3pop_lib.R")
source("~/Packages/s_lattice/simulation.R")

#Modern GBR, CEU, FIN, TSI
#Ancient WHG, ENeo, Yamnaya

########################################################################
## --args: x v8 2

########################################################################

source("~/selection/code/analysis/setup_populations_reads.R")

########################################################################

if(version=="v6"){
    lambda=1.21
    sig <- 10^-6.79                                      #genome-wide significance level
}else if(version=="v8"){
    num.res.tag <- as.numeric(results.tag)
    lambda=c(1.3, 1.3, 1.3)[num.res.tag]
    sig <- 10^-7.30                                      #genome-wide significance level
}

########################################################################

rnds <- seq(0,1,length.out=11)                               #proportions of randomness
Nsims <- 1000                                           #Number of replicates per test
Npowersims <- 1000                                           #Number of replicates per test
gens <- 100
s <- 0.02
Ne <- 14000

########################################################################

selpops <- c( "CEU", "GBR", "IBS", "TSI")
degf <- dim(A)[2]

## Prepare counts and totals
counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
data <- counts[,1:5]
counts <- data.matrix(counts[,6:NCOL(counts)])
totals <- data.matrix(totals[,6:NCOL(totals)])

## setup for read data. 
include.read.samples <- read.samples(indfile, include.reads)
empty.data <- make.empty.data(pops)

#Sample uniformly per chromosome. 
counts.per.chr <- table(sample(data$CHR[data$CHR<=22], Nsims))

#Select sites with f0 < 0.1
#Now these are lists of frequencies. 
tf <- list()
i=1
for(chr in 1:22){
    cat(paste0("chr", chr))
    if(!(as.character(chr) %in% names(counts.per.chr))){
        next
    }

    reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
    k=1
    while(k <= counts.per.chr[as.character(chr)]){
        cat(paste0("\rchr", chr, " ", k, "/", counts.per.chr[as.character(chr)]))
        
        inc <- data$CHR==chr
        try <- sample(sum(inc), 1)
        snp <- data[inc,][try,"ID"]
        this.reads <- reads[reads[,1]==snp,]
        freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts, this.reads, counts[inc,][try,], totals[inc,][try,], empty.data)

        monomorphic <- all(counts[inc,][try,monocheck]==0)|all(counts[inc,][try,monocheck]==totals[inc,][try,monocheck])
        if(monomorphic){next}

        ## model fit gives us frequency of ref allele. 
        ## fr <- 1-fit.unconstrained.model.reads(freq.data, error.prob=error.prob)$par
        ## if(mean(fr)>0.1){next}
        tf[[i]] <- freq.data
        
        i <- i+1
        k <- k+1
    }
}

lambda.all <- rep(0, length(rnds))
#Estimate lambda as a function of r
for( rndi in 1:length(rnds)){
    rnd <- rnds[rndi]
    cat(paste0(rnd, "\n"))
    this.res <- matrix(0, nrow=Nsims, ncol=2)
    
    for(i in 1:Nsims){
        Atmp <- A
        for(j in 1:NCOL(A)){
            kg <- rexp(3)
            rndp <- kg/sum(kg)
            Atmp[,j] <- rnd*rndp+(1-rnd)*A[,j]
        }

        this.res[i,] <- test.3pop.reads(tf[[i]], Atmp, error.prob=error.prob)
    }
    lambda.all[rndi] <- median(this.res[this.res[,1]>0,1])/qchisq(0.5, df=degf)
}

pdf("~/selection/analysis/power/reads_robust_lambda.pdf")
plot(rnds, lambda.all, col="#377EBA", type="b", pch=16, bty="n", lwd=2, xlab="Random proportion", ylab="Genomic inflation factor", ylim=c(1.2, 1.4))
dev.off()

###########################################################################################
#Now test power.
gens <- 200
s <- 0.01
Ne <- 6000
counts.per.chr <- table(sample(data$CHR[data$CHR<=22], Npowersims))

tf <- list()
i=1
for(chr in 1:22){
    cat(paste0("chr", chr))
    reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts"), as.is=TRUE, header=FALSE)
    k=1
    while(k <= counts.per.chr[chr]){
        cat(paste0("\rchr", chr, " ", k, "/", counts.per.chr[chr]))
        inc <- data$CHR==chr
        try <- sample(sum(inc), 1)
        snp <- data[inc,][try,"ID"]
        this.reads <- reads[reads[,1]==snp,]
        freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts, this.reads, counts[inc,][try,], totals[inc,][try,], empty.data)

        monomorphic <- all(counts[inc,][try,monocheck]==0)|all(counts[inc,][try,monocheck]==totals[inc,][try,monocheck])
        if(monomorphic){next}

        ## model fit gives us frequency of ref allele. 
        fr <- 1-fit.unconstrained.model.reads(freq.data, error.prob=error.prob)$par
        if(mean(fr)>0.2){next}
        tf[[i]] <- freq.data
        
        i <- i+1
        k <- k+1
    }
}

one.power <- c(0, length(rnds))
for( rndi in 1:length(rnds)){
    rnd <- rnds[rndi]
    cat(paste0(rnd, "\n"))

    results.one <- 0
    this.tf <- tf
    for(i in 1:Npowersims){
        pop <- sample(selpops, 1)
        this.fr <- pmax(0.01, this.tf[[i]][[pop]][["counts"]][2]/sum(this.tf[[i]][[pop]][["counts"]]))
        traj <- simulate.wright.fisher(Ne, gens, this.fr, s)
        new.fr <- rev(traj)[1]
        this.tot <- sum(this.tf[[i]][[pop]][["counts"]])
        this.alt <-  rbinom(1, this.tot, new.fr)
        this.tf[[i]][[pop]][["counts"]] <- c(this.tot-this.alt, this.alt)
        
        ## Randomise A
        Atmp <- A
        for(j in 1:NCOL(A)){
            kg <- rexp(3)
            rndp <- kg/sum(kg)
            Atmp[,j] <- rnd*rndp+(1-rnd)*A[,j]
        }
        
        test <- test.3pop.reads(this.tf[[i]], Atmp, error.prob=error.prob)
        corrected.p <- pchisq(test[1]/lambda.all[rndi], df=degf, lower.tail=F)
        if(corrected.p<sig){results.one <- results.one+1} 
    }
    one.power[rndi] <- results.one/Npowersims

}
dev.new()
pdf("~/selection/analysis/power/reads_robust_power.pdf")
plot(rnds, one.power, col="#E41A1C", type="b", pch=16, bty="n", lty=1, lwd=2, xlab="Random proportion", ylab="Power")
dev.off()

pdf("~/selection/analysis/power/reads_robust_power_lambda.pdf")
par(mar=c(5,4,4,4))
plot(rnds, lambda.all, col="#377EBA", type="b", pch=16, bty="n", lwd=2, xlab="Random proportion", ylab="Genomic inflation factor", yaxt="n", xaxt="n", ylim=c(1.2,1.4))
axis(1, lwd=2)
axis(2, col="#377EBA", lwd=2)
par(new=TRUE)
plot(rnds, one.power, col="#CC5500", type="b", pch=16, bty="n", lwd=2, axes=FALSE, xlab="", ylab="")
axis(4, col="#CC5500", lwd=2)
mtext("Power", 4, line=3)
dev.off()
