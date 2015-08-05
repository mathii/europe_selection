source("~/selection/code/lib/3pop_lib.R")
source("~/Packages/s_lattice/simulation.R")

########################################################################
## --args: x v8 2 Proportion 12345
## This tests how lambda and power chage under different scenarios,
## Either when the admixture matrix A is misspecified, or if there is
## admixture from YRI into one of the populations. It first selects a
## bunch of loci, then runs through twice - once with neutral loci
## to do the simulation and estimate lambda, and then once more to
## estimate power. 
########################################################################

source("~/selection/code/analysis/setup_populations_reads.R")

########################################################################

if(version=="v6"){
    sig <- 10^-6.79                                      #genome-wide significance level
}else if(version=="v8"){
    num.res.tag <- as.numeric(results.tag)
    sig <- 10^-7.30                                      #genome-wide significance level
}

########################################################################
cA <- commandArgs(TRUE)
which.test <-cA[4]
if(!(which.test %in% c("Proportion", "Admixture"))){
  stop("Test (arg 4) must be Proportions or Admixture")
}
seed <- 12345
if(length(cA>4)){
  seed <- cA[5]
  set.seed(seed)
}
rnds <- seq(0,1,length.out=11)                               #proportions of randomness
if(which.test=="Admixture"){
  rnds <- seq(0,0.1,length.out=11)      #Proportions of admixture
}

########################################################################

randomize.A <- function(A, rnd){
  Atmp <- A
  for(j in 1:NCOL(A)){
    kg <- rexp(3)
    rndp <- kg/sum(kg)
    Atmp[,j] <- rnd*rndp+(1-rnd)*A[,j]
  }
  return(Atmp)
}

########################################################################
add.YRI.admixture <- function(this.tf, rnd, selpops, counts, totals){
  for(i in 1:length(this.tf)){
    snp=names(this.tf)[i]
    adm <- sample(selpops,1)
    ttfc <- this.tf[[i]][[adm]][["counts"]]
    nt <- sum(ttfc)
    sub.ref.alt <- round((1-rnd)* ttfc)
    to.fill <- nt-sum(sub.ref.alt)  #total number to pick from YRI
    fill.c <- round(to.fill*counts[which(data$ID==snp),"YRI"]/totals[which(data$ID==snp),"YRI"])
    this.tf[[i]][[adm]][["counts"]] <- sub.ref.alt+c(fill.c, to.fill-fill.c)
  }
  return(this.tf)
}

########################################################################

Npowersims <- Nlambdasims <- 1000                            #Number of replicates per test
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
rownames(counts) <- data$ID
rownames(totals) <- data$ID

## setup for read data. 
include.read.samples <- read.samples(indfile, include.reads)

max.freq <- 1                           #Actually max.freq. 
tf.file <- paste0("~/selection/analysis/", version,"/power/tf_",which.test,"_seed_", seed, "_maxf_", max.freq, "_N_", Nlambdasims,".obj")
if(file.exists(tf.file)){
  load(tf.file)
  cat(paste0("Loading seed file ", tf.file, "\n"))
}else{
  tf <- sample.data(data, Nlambdasims, read.root, pops, include.reads, include.read.samples, include.counts, counts, totals, max.freq, monocheck)
  save(tf, file=tf.file)
}

lambda.all <- rep(0, length(rnds))
#Estimate lambda as a function of r
for( rndi in 1:length(rnds)){
  rnd <- rnds[rndi]
  cat(paste0(rnd, "\n"))
  this.res <- matrix(0, nrow=Nlambdasims, ncol=2)
  
  this.tf <- tf
  if(which.test=="Admixture"){
    this.tf <- add.YRI.admixture(tf, rnd, selpops, counts, totals)
  }
  for(i in 1:Nlambdasims){
    Atmp <- A
    if(which.test=="Proportion"){
      Atmp <- randomize.A(A, rnd)
    }
    this.res[i,] <- test.3pop.reads(this.tf[[i]], Atmp, error.prob=error.prob)
  }
  lambda.all[rndi] <- median(this.res[this.res[,1]>0,1])/qchisq(0.5, df=degf)
}

## pdf(paste0("~/selection/analysis/",version,"/power/reads_robust_", which.test,"_seed_", seed,"_lambda.pdf"))
## plot(rnds, lambda.all, col="#377EBA", type="b", pch=16, bty="n", lwd=2, xlab="Random proportion", ylab="Genomic inflation factor", ylim=c(1.2, 1.4))
## dev.off()

###########################################################################################
#Now test power.
## Here we're restricting to things with a MAF < 0.1

max.freq <- 0.1
tf.file <- paste0("~/selection/analysis/", version,"/power/tf_",which.test,"_seed_", seed, "_maxf_", max.freq, "_N_", Nlambdasims,".obj")
if(file.exists(tf.file)){
  load(tf.file)
  cat(paste0("Loading seed file ", tf.file, "\n"))
}else{
  tf <- sample.data(data, Nlambdasims, read.root, pops, include.reads, include.read.samples, include.counts, counts, totals, max.freq, monocheck)
  save(tf, file=tf.file)
}

all.power <- rep(0, length(rnds))
for( rndi in 1:length(rnds)){
    rnd <- rnds[rndi]
    cat(paste0(rnd, "\n"))

    results.all <- 0
    this.tf <- tf
    if(which.test=="Admixture"){
      this.tf <- add.YRI.admixture(tf, rnd, selpops, counts, totals)
    }
    for(i in 1:Npowersims){
        for(pop in selpops){
          this.fr <- pmax(0.01, this.tf[[i]][[pop]][["counts"]][2]/sum(this.tf[[i]][[pop]][["counts"]]))
          traj <- simulate.wright.fisher(Ne, gens, this.fr, s)
          new.fr <- rev(traj)[1]
          this.tot <- sum(this.tf[[i]][[pop]][["counts"]])
          this.alt <-  rbinom(1, this.tot, new.fr)
          this.tf[[i]][[pop]][["counts"]] <- c(this.tot-this.alt, this.alt)
        }
        
        ## Randomise A
        Atmp <- A
        if(which.test=="Proportion"){
          Atmp <- randomize.A(A, rnd)
        }

        test <- test.3pop.reads(this.tf[[i]], Atmp, error.prob=error.prob)
        corrected.p <- pchisq(test[1]/lambda.all[rndi], df=degf, lower.tail=F)
        if(corrected.p<sig){results.all <- results.all+1} 
    }
    all.power[rndi] <- results.all/Npowersims

}

results <- data.frame(random=rnds, lambda=lambda.all, power=all.power)
write.table(results, paste0("~/selection/analysis/", version,"/power/reads_robust_power_", which.test,"_seed_", seed,"_lambda.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## pdf(paste0("~/selection/analysis/", version,"/power/reads_robust_power_", which.test,"_seed_", seed,"_lambda.pdf"))
## par(mar=c(5,4,4,4))
## plot(rnds, lambda.all, col="#377EBA", type="b", pch=16, bty="n", lwd=2, xlab="Random proportion", ylab="Genomic inflation factor", yaxt="n", xaxt="n", ylim=c(1.2,1.4))
## axis(1, lwd=2)
## axis(2, col="#377EBA", lwd=2)
## par(new=TRUE)
## plot(rnds, all.power, col="#CC5500", type="b", pch=16, bty="n", lwd=2, axes=FALSE, xlab="", ylab="")
## axis(4, col="#CC5500", lwd=2)
## mtext("Power", 4, line=3)
## dev.off()
