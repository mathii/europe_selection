#Test whether the modern population frequencies can be modelled as a mixture of the
#Three ancestral populations.
library(reshape2)
source("~/selection/code/lib/3pop_lib.R")
dir <- getwd()
setwd("~/Packages/s_lattice/")
source("include.R")
setwd(dir)

#Modern GBR, CEU, FIN, TSI
#Ancient WHG, ENeo, Yamnaya

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

gss <- c(50, 100, 200)                               #Generations of selection
## ss <- 10^(seq(log10(0.002), log10(0.1), length.out=10)) #Selection coefficient
ss <- c(0.002, 0.003, 0.005, 0.08, 0.01, 0.02, 0.03, 0.05, 0.08, 0.1)
Ne <- 6000                                          #2 Population size
N <- 1000                                           #Number of replicates

########################################################################

selpops <- c( "CEU", "GBR", "IBS", "TSI")
#Check if the SNP is monomorphic in these populations. 
A <- matrix(c(0.164, 0.366, 0.470, 0.213, 0.337, 0.450, 0, 0.773, 0.227, 0, 0.712, 0.288),3, 4) 
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
counts.per.chr <- table(sample(data$CHR[data$CHR<=22], N))

#Select sites with f0 < 0.2
#Now these are lists of frequencies. 
tf <- list()
i=1
for(chr in 1:22){
    cat(paste0("\rchr", chr))
    if(!(as.character(chr) %in% names(counts.per.chr))){
        next
    }
    
    reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
    k=1

    while(k <= counts.per.chr[as.character(chr)]){
        cat(paste0("\rchr", chr, " ", k, "/", counts.per.chr[chr]))
        inc <- data$CHR==chr
        try <- sample(sum(inc), 1)
        snp <- data[inc,][try,"ID"]
        this.reads <- reads[reads[,1]==snp,]
        freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts, this.reads, counts[inc,][try,], totals[inc,][try,], empty.data)


        ## model fit gives us frequency of ref allele. 
        fr <- 1-fit.unconstrained.model.reads(freq.data, error.prob=error.prob)$par
        if(mean(fr)>0.2){next}

        ## Don't want it to me monomporphic. 
        monomorphic <- all(counts[inc,][try,monocheck]==0)|all(counts[inc,][try,monocheck]==totals[inc,][try,monocheck])
        if(monomorphic){next}

        tf[[i]] <- freq.data
        
        i <- i+1
        k <- k+1
    }
}

results.all <- matrix(0,nrow=length(ss), ncol=length(gss))
results.one <- matrix(0,nrow=length(ss), ncol=length(gss))
#Add selection (to all pops)
for( gsi in 1:length(gss)){
    gs <- gss[gsi]
    for(si in 1:length(ss)){
        cat(paste0(gss[gsi]," ",  ss[si], "\n"))

        #All pops
        this.tf <- tf
        for(i in 1:N){
            for(pop in selpops){
                this.fr <- pmax(0.01, this.tf[[i]][[pop]][["counts"]][2]/this.tf[[i]][[pop]][["counts"]][1])
                traj <- simulate.wright.fisher(Ne, gss[gsi], this.fr, ss[si])
                new.fr <- rev(traj)[1]
                this.tot <- sum(this.tf[[i]][[pop]][["counts"]])
                this.alt <-  rbinom(1, this.tot, new.fr)
                this.tf[[i]][[pop]][["counts"]] <- c(this.tot-this.alt, this.alt)
            }
            test <- test.3pop.reads(this.tf[[i]], A, error.prob=error.prob)
            test[2] <- pchisq(test[1]/lambda, df=degf, lower.tail=FALSE)
            if(test[2]<sig){results.all[si,gsi] <- results.all[si,gsi]+1} 
        }

        #One pop
        this.tf <- tf
        for(i in 1:N){
            pop <- sample(selpops, 1)
            this.fr <- pmax(0.01, this.tf[[i]][[pop]][["counts"]][2]/this.tf[[i]][[pop]][["counts"]][1])
            traj <- simulate.wright.fisher(Ne, gss[gsi], this.fr, ss[si])
            new.fr <- rev(traj)[1]
            this.tot <- sum(this.tf[[i]][[pop]][["counts"]])
            this.alt <-  rbinom(1, this.tot, new.fr)
            this.tf[[i]][[pop]][["counts"]] <- c(this.tot-this.alt, this.alt)
            test <- test.3pop.reads(this.tf[[i]], A, error.prob=error.prob)
            test[2] <- pchisq(test[1]/lambda, df=degf, lower.tail=FALSE)
            if(test[2]<sig){results.one[si,gsi] <- results.one[si,gsi]+1} 
        }

    }
}
results.all <- results.all/N
results.one <- results.one/N


pdf(paste0("~/selection/analysis/",version,"/power/read_power", results.tag,".pdf"))
plot(ss, results.all[,1], col="#377EBA", type="b", pch=16, bty="n", lty=2, ylim=c(0,1), log="x", xlab="Selection coefficient", ylab="power")
lines(ss, results.all[,2], col="#E41A1C", type="b", pch=16, lty=2)
lines(ss, results.all[,3], col="#4DAF4A", type="b", pch=16, lty=2)
lines(ss, results.one[,1], col="#377EBA", type="b", pch=1, lty=3)
lines(ss, results.one[,2], col="#E41A1C", type="b", pch=1, lty=3)
lines(ss, results.one[,3], col="#4DAF4A", type="b", pch=1, lty=3)
legend("topleft", c("Selected in all populations", "Selected in one population", "50 generations of selection", "100 generations of selection", "200 generations of selection"), col=c("black", "black", "#377EBA", "#E41A1C", "#4DAF4A"), pch=16, lty=c(2,3,1,1,1), bty="n")
dev.off()

colnames(results.all) <- colnames(results.one) <- gss
rownames(results.all) <- rownames(results.one) <- ss
m1 <- cbind(pop="All", melt(results.all))
m2 <- cbind(pop="One", melt(results.one))
mres <- rbind(m1, m2)
if(results.tag!=""){
    mres <- cbind(results.tag, mres)
    names(mres)[1]<-"Tag"
}
mres <- cbind(mres, round(sum(rowMeans(sapply(tf,  effective.data.size)[1:3,])), 1))
names(mres)[NCOL(mres)] <- "Eff.N"
write.table(mres,paste0("~/selection/analysis/",version,"/power/read_power", results.tag,".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
