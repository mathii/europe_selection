#Functions to implement the 3 (or n) population test.
## require(numDeriv)

EPSILON.3pop <- 1e-4 #Minumum allele frequency

#########################################################
#
# Likelihood of allele counts N.A/N given population
# frequencies p in source and Ap in derived pops
#
#########################################################

constrained.likelihood <- function(p, A, N, N.A){
    p <- c(p, c(p %*% A))
    p <- pmin(pmax(EPSILON.3pop, p), 1-EPSILON.3pop)
    return(sum(dbinom(N.A, N, p, log=TRUE)))
}

#########################################################
#
# Fit the allele frequencies in the source populations. 
#                                        
#########################################################


fit.constrained.model <- function(N, N.A, A){
    p.anc.init <- rep(0.5, dim(A)[1])
    opt <- optim(p.anc.init, constrained.likelihood, A=A, N=N, N.A=N.A, control=list(fnscale=-1), lower=0, upper=1, method="L-BFGS")
    return(opt)
}

#########################################################
#
# Likelihood ratio test, given observations N.A/N and
# matrix A mapping source to derived populations.
# N and N.A are length K and A is a matrix with
# sum(dim(A))==K
#
#########################################################

test.3pop <- function(N, N.A, A){
    degf <- dim(A)[2]

    if(length(N)!=length(N.A)){stop("N and N.A different lengths")}
    if(sum(dim(A))!=length(N)){stop("Matrix A not compatible with observations")}

    #Alternative
    #MLE is just the observed proportion
    #If no observations then it doesn't matter what p is
    p.hat <- ifelse(N==0, 0.5, N.A/N)  
    p.hat <- pmin(pmax(EPSILON.3pop, p.hat), 1-EPSILON.3pop)
    l1 <- sum(dbinom(N.A, N, p.hat, log=TRUE))

    #Null
    opt <- fit.constrained.model(N, N.A, A)
    p.est <- c(opt$par, opt$par %*% A) 
    p.est <- pmin(pmax(EPSILON.3pop, p.est), 1-EPSILON.3pop)
    l0 <- sum(dbinom(N.A, N, p.est, log=TRUE))

    stat <- 2*(l1-l0)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}

#########################################################
#
# Likelihood ratio test, given observations N.A/N, test
# whether all the observations can come from the
# same population. 
#
#########################################################

test.diff <- function(N, N.A){
    degf <- length(N)-1

    if(length(N)!=length(N.A)){stop("N and N.A different lengths")}

    #Alternative
    #MLE is just the observed proportion
    #If no observations then it doesn't matter what p is
    p.hat <- ifelse(N==0, 0.5, N.A/N)  
    p.hat <- pmin(pmax(EPSILON.3pop, p.hat), 1-EPSILON.3pop)
    l1 <- sum(dbinom(N.A, N, p.hat, log=TRUE))

    #Null
    p.est <- rep(sum(N.A)/sum(N), length(N))
    p.est <- pmin(pmax(EPSILON.3pop, p.est), 1-EPSILON.3pop)
    l0 <- sum(dbinom(N.A, N, p.est, log=TRUE))

    stat <- 2*(l1-l0)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}

#########################################################
#
# Likelihood for the combined read and count based
# model.
# p: Reference allele frequencies (length N)
# data: list of N populations, each entry is a list with
# three elements, "reads", "probabilities" and "counts". "reads" is itself
# a list with "ref" and "alt" vectors of ref and alt
# read counts for each individual and "counts" is a vector
# of length 2 for giving ref/alt counts. Probabilities is a list of,
# genotype probabilities (00, 01, and 11), named by samples. 
# This returns the log likelihood summed over all the
# populations. het.p is the probability of seeing the reference
# read for heterozygotes (0.5 by default)
#########################################################

likelihood.reads <- function(freq, data, error.prob=0.001, het.p=0.5){
    ll <- 0
    i.pop <- 1
    for(pop in names(data)){
        N.read.ind <- length(data[[pop]][["reads"]][["ref"]]) #Number of individuals with read information
        p <- freq[i.pop]
        ## If the frequency is NA, and there's no data; likelihood change is zero
        if(is.na(p) & (N.read.ind==0) & (sum(data[[pop]][["counts"]])==0)){
          i.pop <- i.pop+1              #Have to increment though.
          next
        }
        if(N.read.ind){
            #Add reads
            for(i in 1:N.read.ind){
                ref <- data[[pop]][["reads"]][["ref"]][[i]]
                alt <- data[[pop]][["reads"]][["alt"]][[i]]
                ll <- ll+log(dbinom(alt, ref+alt, error.prob)*p*p + dbinom(ref, ref+alt, error.prob)*(1-p)*(1-p) + dbinom(ref, ref+alt, het.p)*2*p*(1-p))
            }
        }
        #Add probabilities
        if(length(data[[pop]][["probabilities"]])>0){
            for(ps in data[[pop]][["probabilities"]]){
                #Is this right?
                pse <- pmin(pmax(ps, error.prob), 1-error.prob)
                ll <- ll+log(pse[1]*p*p + pse[2]*2*p*(1-p) + pse[3]*(1-p)*(1-p))
            }
        }
        #Add counts, can be zero
        if("counts" %in% names(data[[pop]])){
            cnts <- data[[pop]][["counts"]]
            ll <- ll+dbinom(cnts[1], cnts[1]+cnts[2], p, log=TRUE)
        }
        i.pop <- i.pop+1    
    }
    return(ll)
}

#########################################################
#
# Fit the constrained model
# 
#########################################################

fit.fixed.model.reads <- function(data, error.prob=0){
    p.anc.init <- 0.5
    opt <- optimize(fixed.likelihood.reads, lower=0, upper=1, data=data, error.prob=error.prob, maximum=TRUE)
    par=list(par=rep(opt$maximum, length(data)), value=opt$objective)
    return(par)
}

#########################################################
#
# Likelihood for read based model with all frequencies
# the same
# 
#########################################################

fixed.likelihood.reads <- function(p, data, error.prob=0){
    p <- rep(p, length(data))
    p <- pmin(pmax(EPSILON.3pop, p), 1-EPSILON.3pop)
    ll <- likelihood.reads(p, data, error.prob=error.prob)
    return(ll)
}

#########################################################
#
# Constrained likelihood for read based model
# 
#########################################################

constrained.likelihood.reads <- function(p, A, data, error.prob=0){
    p <- c(p, c(p %*% A))
    p <- pmin(pmax(EPSILON.3pop, p), 1-EPSILON.3pop)
    ll <- likelihood.reads(p, data, error.prob=error.prob)
    return(ll)
}

#########################################################
#
# Fit the constrained model
# 
#########################################################

fit.constrained.model.reads <- function(data, A, error.prob=0){
    p.anc.init <- rep(0.5, dim(A)[1])
    opt <- optim(p.anc.init, constrained.likelihood.reads, A=A, data=data, error.prob=error.prob, control=list(fnscale=-1), lower=0, upper=1, method="L-BFGS")
    return(opt)
}

#########################################################
#
# Fit the unconstrained model
# Since the frequencies are independent, we can fit them
# separately. 
# 
#########################################################

fit.unconstrained.model.reads <- function(data, error.prob=0.001){
    p <- rep(NA, length(data))
    for(i in 1:length(data)){
        subdata=list(data[[i]])
        names(subdata) <- names(data)[i]
        if(is.null(subdata[[1]]$reads$ref)&(sum(subdata[[1]]$counts)==0)){next} #If there's no data at all. 
        opt <- optimize(likelihood.reads, data=subdata, error.prob=error.prob, maximum=TRUE, lower=EPSILON.3pop, upper=1-EPSILON.3pop)
        p[i] <- opt$maximum
    }
    ll <- likelihood.reads(p, data, error.prob=error.prob)
    return(list(par=p, value=ll))
}

#########################################################
#
# Helper function. 
#
#########################################################

ci.optfun <- function(f, data, error.prob, target){
    return(target-likelihood.reads(f, data, error.prob=error.prob))
}

#########################################################
#
# Upper confidence intervals for unconstrained model 
#
#########################################################

ci.unconstrained.model.reads <- function(data, error.prob=0.001, alpha=0.95){
    uci <- rep(NA, length(data))
    lci <- rep(NA, length(data))
    f <- rep(NA, length(data))
    ll <- rep(NA, length(data))
    ll.diff <- qchisq(alpha, df=1)/2
    for(i in 1:length(data)){
        subdata=list(data[[i]])
        names(subdata) <- names(data)[i]

        #If there's no data at all. 
        if(is.null(subdata[[1]]$reads$ref)&(sum(subdata[[1]]$counts)==0)){
          ll[i] <- 0
          next
        }
        opt <- optimize(likelihood.reads, data=subdata, error.prob=error.prob, maximum=TRUE, lower=EPSILON.3pop, upper=1-EPSILON.3pop)
        f[i] <- opt$maximum
        ll[i] <- likelihood.reads(f[i], subdata, error.prob=error.prob)
                    
        #Upper
        if(ll[i]-likelihood.reads(1-EPSILON.3pop, subdata, error.prob=error.prob)<ll.diff){
            uci[i] <- 1
        }else{
            uci[i] <- uniroot(ci.optfun, interval=c(f[i], 1-EPSILON.3pop), data=subdata, error.prob=error.prob, target=ll[i]-ll.diff)$root
        }

        #Lower
        if(ll[i]-likelihood.reads(EPSILON.3pop, subdata, error.prob=error.prob)<ll.diff){
            lci[i] <- 0
        }else{
            lci[i] <- uniroot(ci.optfun, interval=c(EPSILON.3pop, f[i]), data=subdata, error.prob=error.prob, target=ll[i]-ll.diff)$root
        }
    }

    return(list(p=f, ll=ll, uci=uci, lci=lci))
}

#########################################################
#
# Effective size of data - roughly how many hard                                       
# calls is this data equivalent to? 
#########################################################

effective.data.size <- function(data){
    data.size <- rep(0, length(data))
    names(data.size) <- names(data)
    for(pop in names(data)){
        total <- 0
        N.read.ind <- length(data[[pop]][["reads"]][["ref"]]) #Number of individuals with read information
        if(N.read.ind){
            #Add reads
            for(i in 1:N.read.ind){
                ref <- data[[pop]][["reads"]][["ref"]][[i]]
                alt <- data[[pop]][["reads"]][["alt"]][[i]]
                total <- total+2-0.5^(ref+alt-1)
            }
        }
        #Add counts, can be zero
        if("counts" %in% names(data[[pop]])){
            cnts <- data[[pop]][["counts"]]
            total <- total+cnts[1]+cnts[2]
        }
        data.size[pop] <- total
    }
    return(data.size)
}

#########################################################
#
# General 3 population test using read information                                        
# 
#########################################################

test.3pop.reads <- function(data, A, error.prob=error.prob){
    degf <- dim(A)[2]

    if(sum(dim(A))!=length(data)){stop("Matrix A not compatible with observations")}

    uf <- fit.unconstrained.model.reads(data, error.prob=error.prob)
    cf <- fit.constrained.model.reads(data, A, error.prob=error.prob)   
    
    stat <- 2*(uf$value-cf$value)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}

#########################################################
#
# Test for a difference in population means, using read                                        
# 
#########################################################

test.diff.reads <- function(data, error.prob=error.prob){
    degf=length(data)-1
    
    uf <- fit.unconstrained.model.reads(data, error.prob=error.prob)
    cf <- fit.fixed.model.reads(data, error.prob=error.prob)
    
    stat <- 2*(uf$value-cf$value)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}

#########################################################
#
# Estimated frequency and effective sample size for
# one data point. Freq.data of size 1 
#########################################################

effective.size.reads <- function(data, error.prob=0){
    uf <- fit.unconstrained.model.reads(data, error.prob=error.prob)
    d2 <- hessian(likelihood.reads, uf$par, method="Richardson", method.args=list(eps=EPSILON.3pop/2), data, error.prob=error.prob)
    N <- -uf$par*(1-uf$par)*c(d2)
    return(list(p=uf$par, N=N, d2=c(d2)))
}

#########################################################
#
# Empty read/allele count data structure
#
#########################################################

make.empty.data <- function(pops){
    empty.data <- rep( list(list()), length(pops) ) 
    names(empty.data) <- pops
    for(pop in pops){
        empty.data[[pop]] <- list("reads"=list("ref"=NULL, "alt"=NULL, "samples"=NULL),
                                  "counts"=c(0,0),
                                  "probabilities"=list()) #ref and alt counts. 
    }
    return(empty.data)
}

#########################################################
#
# Combine information from read counts and allele counts 
# into the data structure used for the likelihood
# calcluations here
#
# pops: list of populations to include
# include.reads: list, with names in pops, mapping to
# vectors of sub-populations to add together to make
# each population, using read counts.
# include.read samples: which samples are in each sub-population
# include.counds: As include.reads, but with populations
# for which we should add hard counts.
# read.info, counts, totals: standard data structures
#########################################################

make.freq.data <- function(pops, include.reads, include.read.samples, include.counts, read.info, this.counts, this.totals, empty.data, include.prob.samples, prob.data){
    freq.data <- empty.data
    for(pop in pops){
        if(pop %in% names(include.reads)){
            for(sample in include.read.samples[[pop]]){
                ref.alt <- read.info[read.info[,2]==sample,3:4]

                if(!(sample %in% read.info[,2])){stop(paste0(sample, "not in read info"))}
                
                if(sum(ref.alt)>0){
                    freq.data[[pop]][["reads"]][["ref"]] <- c(freq.data[[pop]][["reads"]][["ref"]],ref.alt[[1]])
                    freq.data[[pop]][["reads"]][["alt"]] <- c(freq.data[[pop]][["reads"]][["alt"]],ref.alt[[2]])
                    freq.data[[pop]][["reads"]][["samples"]] <- c(freq.data[[pop]][["reads"]][["samples"]],sample)

                }
            }
        }

        if(pop %in% names(include.counts)){
            for(subpop in include.counts[[pop]]){
                ref.alt <- c(this.counts[subpop][[1]], this.totals[subpop][[1]]-this.counts[subpop][[1]])
                names(ref.alt) <- NULL
                freq.data[[pop]][["counts"]] <- freq.data[[pop]][["counts"]]+ref.alt
            }
        }
    }
    return(freq.data)
}

#########################################################
#
# Combine information from read counts and genotype probabilities 
# into the data structure used for the likelihood
# calcluations here
#
# pops: list of populations to include
# include.reads: list, with names in pops, mapping to
# vectors of sub-populations to add together to make
# each population, using read counts.
# include.read samples: which samples are in each sub-population
# include.counds: As include.reads, but with populations
# for which we should add hard counts.
# read.info, counts, totals: standard data structures
#########################################################

make.prob.freq.data <- function(pops, include.probs, include.prob.samples, include.counts, prob.info, this.counts, this.totals, empty.data){
    freq.data <- empty.data
    for(pop in pops){
        if(pop %in% names(include.probs)){
            for(sample in include.prob.samples[[pop]]){
                this.str <- as.character(prob.info[sample])
                this.gp <- as.numeric(strsplit( strsplit(this.str, ":", fixed=TRUE)[[1]][3], ",")[[1]])
                freq.data[[pop]][["probabilities"]][[sample]] <- this.gp
            }
        }

        if(pop %in% names(include.counts)){
            for(subpop in include.counts[[pop]]){
                ref.alt <- c(this.counts[subpop][[1]], this.totals[subpop][[1]]-this.counts[subpop][[1]])
                names(ref.alt) <- NULL
                freq.data[[pop]][["counts"]] <- freq.data[[pop]][["counts"]]+ref.alt
            }
        }
    }
    return(freq.data)
}


#########################################################
#
## get list of samples in each population of reads
## 
#########################################################

read.samples <- function(indfile, include.reads, exclude=c()){
    ind <- read.table(indfile, as.is=TRUE, header=FALSE)
    include.read.samples <- lapply(include.reads, function(x){NULL})
    for(pop in names(include.reads)){
        for(subpop in include.reads[[pop]]){
            to.include <- ind[ind[,3]==subpop,1]
            to.include <- to.include[!(to.include %in% exclude)]
            include.read.samples[[pop]] <- c(include.read.samples[[pop]], to.include)
        }
    }
    return(include.read.samples)
}

#########################################################
## 
## Sample some number of reads, uniformly distributed
## over the autosomes, for power simulations. 
## 
#########################################################

sample.data <- function(data, N, read.root, pops, include.reads, include.read.samples, include.counts, counts, totals, max.freq, monocheck){
  counts.per.chr <- table(sample(data$CHR[data$CHR<=22], N))
  empty.data <- make.empty.data(pops)

  tf <- list()
  i=1
  for(chr in 1:22){
    cat(paste0("\rchr", chr))
    if(!(as.character(chr) %in% names(counts.per.chr))){
      next
    }
    
    reads <- read.table(paste0(read.root, ".chr", chr, ".readcounts.gz"), as.is=TRUE, header=FALSE)
    k=1
    inc <- data$CHR==chr
    total.inc <- sum(inc)
    inc.counts <- counts[inc,]
    inc.totals <- totals[inc,]
    inc.data <- data[inc,]
    while(k <= counts.per.chr[as.character(chr)]){
        cat(paste0("\rchr", chr, " ", k, "/", counts.per.chr[chr]))
        try <- sample(total.inc, 1)
        snp <- inc.data[try,"ID"]
        this.reads <- reads[reads[,1]==snp,]
        freq.data <- make.freq.data(pops, include.reads, include.read.samples, include.counts, this.reads,  inc.counts[try,], inc.totals[try,], empty.data)

        ## model fit gives us frequency of ref allele. 
        fr <- 1-fit.unconstrained.model.reads(freq.data, error.prob=error.prob)$par
        if(mean(fr, na.rm=TRUE)>max.freq){next}

        ## Don't want it to me monomporphic. 
        monomorphic <- all(inc.counts[try,monocheck]==0)|all(inc.counts[try,monocheck]==inc.totals[try,monocheck])
        if(monomorphic){next}

        tf[[i]] <- freq.data
        names(tf)[i] <- snp
        
        i <- i+1
        k <- k+1
      }
  }
  return(tf)
}
